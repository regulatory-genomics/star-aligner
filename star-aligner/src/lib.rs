use anyhow::{bail, Result};
use noodles::{
    fastq,
    sam::{
        self,
        header::record::value::{map::ReferenceSequence, Map},
        Header,
    },
};
use rayon::prelude::*;
use std::{
    ffi::{c_char, c_int, CStr, CString},
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    sync::Arc,
    vec,
};

/// Struct representing STAR aligner options.
///
/// `StarOpts` stores settings to control the behavior of the STAR aligner.
#[derive(Clone)]
pub struct StarOpts {
    genome_dir: PathBuf, // Path to the STAR reference genome
    num_threads: u16,    // Number of threads to use
    sjdb_overhang: u16,  // Overhang for splice junctions
    out_filter_score_min: u16, // Alignment will be output only if its score is higher than or equal to this value.
                               // To output higher quality reads, you can set this value closer to read length.
                               // Default: 0.
}

impl StarOpts {
    /// Create a new `StarOpts` instance.
    ///
    /// # Arguments
    /// - `genome_dir` - Path to the STAR genome directory.
    ///
    /// # Returns
    /// - A `StarOpts` object with default values for threads, splice junction
    ///   overhang, and minimum filter score.
    pub fn new<P: AsRef<Path>>(genome_dir: P) -> Self {
        Self {
            genome_dir: genome_dir.as_ref().to_path_buf(),
            num_threads: 8,
            sjdb_overhang: 100,
            out_filter_score_min: 0,
        }
    }

    /// Converts options to a vector of C-style strings for FFI compatibility.
    ///
    /// # Returns
    /// - A vector of pointers to C-style strings representing each option.
    fn to_c_args(&self) -> Vec<*const i8> {
        let num_threads = self.num_threads.to_string();
        let sjdb_overhang = self.sjdb_overhang.to_string();
        let out_filter_score_min = self.out_filter_score_min.to_string();
        let args = vec![
            "STAR",
            "--genomeDir",
            self.genome_dir.to_str().unwrap(),
            "--runThreadN",
            num_threads.as_str(),
            "--sjdbOverhang",
            sjdb_overhang.as_str(),
            "--outFilterScoreMin",
            out_filter_score_min.as_str(),
            "--readNameSeparator",
            "space",
            "--outSAMunmapped",
            "Within",
            "KeepPairs",
            "--outSAMtype",
            "SAM",
            "--outStd",
            "SAM",
            "--outSAMorder",
            "PairedKeepInputOrder",
        ];
        args.into_iter()
            .map(|s| CString::new(s).unwrap().into_raw() as *const i8)
            .collect()
    }
}

/// STAR index representation.
///
/// `StarIndex` stores the reference index and header for STAR alignments.
struct StarIndex {
    index: *const star_sys::StarRef,
    header: Header,
}

unsafe impl Send for StarIndex {}
unsafe impl Sync for StarIndex {}

impl Drop for StarIndex {
    fn drop(&mut self) {
        unsafe {
            star_sys::destroy_ref(self.index);
        }
    }
}

impl StarIndex {
    /// Load a reference index into memory based on the given settings.
    fn load(opts: &StarOpts) -> Result<Self> {
        let header = generate_header(opts.genome_dir.as_path())?;

        let c_args = opts.to_c_args();
        let length = c_args.len() as c_int;

        let index = unsafe { star_sys::init_star_ref(length, c_args.as_ptr()) };
        Ok(StarIndex { index, header })
    }
}

/// A wrapper around the STAR aligner. This is solely used to ensure that the
/// aligner is properly deallocated when it goes out of scope.
pub struct InnerAligner {
    inner: *mut star_sys::Aligner,
    sam_buf: Vec<u8>,
    fq1_buf: Vec<u8>,
    fq2_buf: Vec<u8>,
}

impl InnerAligner {
    fn new(index: &StarIndex) -> Self {
        let inner = unsafe { star_sys::init_aligner_from_ref(index.index) };
        Self {
            inner,
            sam_buf: Vec::new(),
            fq1_buf: Vec::new(),
            fq2_buf: Vec::new(),
        }
    }

    /// Align a single read. Note that the aligner is mutably borrowed.
    /// We did this on purpose to ensure this function cannot be called concurrently.
    pub fn align_read(&mut self, fq: &fastq::Record) -> Result<Vec<sam::Record>> {
        let fq_buf = &mut self.fq1_buf;
        let sam_buf = &mut self.sam_buf;

        // Copy the fastq record into a buffer
        copy_fq_to_buf(fq_buf, fq)?;

        let fq_cstr = CStr::from_bytes_with_nul(fq_buf)?;
        let res: *const c_char = unsafe { star_sys::align_read(self.inner, fq_cstr.as_ptr()) };
        if res.is_null() {
            bail!("STAR returned null alignment");
        }

        let aln_buf = unsafe { CStr::from_ptr(res) }.to_bytes();
        let records: Vec<_> = parse_sam_to_records(aln_buf, sam_buf, fq.name()).collect();

        unsafe {
            libc::free(res as *mut libc::c_void);
        }
        Ok(records)
    }

    /// Align a pair of reads. Note that the aligner is mutably borrowed.
    /// We did this on purpose to ensure this function cannot be called concurrently.
    pub fn align_read_pair(
        &mut self,
        fq1: &fastq::Record,
        fq2: &fastq::Record,
    ) -> Result<(Vec<sam::Record>, Vec<sam::Record>)> {
        let fq_buf1 = &mut self.fq1_buf;
        let fq_buf2 = &mut self.fq2_buf;
        let sam_buf = &mut self.sam_buf;

        // Copy the fastq records into buffers
        copy_fq_to_buf(fq_buf1, fq1)?;
        copy_fq_to_buf(fq_buf2, fq2)?;
        let fq1_str = CStr::from_bytes_with_nul(fq_buf1)?;
        let fq2_str = CStr::from_bytes_with_nul(fq_buf2)?;

        let res: *const c_char =
            unsafe { star_sys::align_read_pair(self.inner, fq1_str.as_ptr(), fq2_str.as_ptr()) };
        if res.is_null() {
            bail!("STAR returned null alignment");
        }

        let aln_buf = unsafe { CStr::from_ptr(res) }.to_bytes();
        let (first_mate, second_mate) = parse_sam_to_records(aln_buf, sam_buf, fq1.name())
            .partition(|rec| rec.flags().unwrap().is_first_segment());

        unsafe {
            libc::free(res as *mut libc::c_void);
        }
        Ok((first_mate, second_mate))
    }
}

impl Drop for InnerAligner {
    fn drop(&mut self) {
        unsafe { star_sys::destroy_aligner(self.inner) };
    }
}

/// Main interface for performing read alignments using STAR.
///
/// `StarAligner` handles alignment of both single and paired-end reads to a reference genome.
pub struct StarAligner {
    index: Arc<StarIndex>,
    opts: StarOpts,
}

impl StarAligner {
    /// Loads the genome index into memory.
    ///
    /// # Arguments
    /// - `opts` - Options controlling alignment parameters.
    ///
    /// # Returns
    /// - A `StarAligner` object initialized with the provided options.
    pub fn new(opts: StarOpts) -> Result<Self> {
        let index = StarIndex::load(&opts)?;
        Ok(Self {
            index: Arc::new(index),
            opts,
        })
    }

    /// Sets the number of threads for alignment.
    ///
    /// # Arguments
    /// - `n` - Number of threads to use.
    ///
    /// # Returns
    /// - `Self` with the updated thread count.
    pub fn with_num_threads(mut self, n: u16) -> Self {
        self.opts.num_threads = n;
        self
    }

    /// Retrieves the header information.
    ///
    /// # Returns
    /// - A reference to the alignment header.
    pub fn get_header(&self) -> &Header {
        &self.index.header
    }

    /// Aligns a batch of single-end reads.
    ///
    /// # Arguments
    /// - `records` - A slice of FASTQ records to align.
    ///
    /// # Returns
    /// - An iterator over vectors of aligned SAM records.
    pub fn align_reads<'a>(
        &'a self,
        records: &'a [fastq::Record],
    ) -> impl ParallelIterator<Item = Vec<sam::Record>> + 'a {
        let chunk_size = self.get_chunk_size(records.len());
        self.init_thread_pool().install(|| {
            records.par_chunks(chunk_size).flat_map_iter(|chunk| {
                let mut aligner = self.clone_inner_aligner();
                chunk
                    .iter()
                    .map(move |fq| aligner.align_read(fq).unwrap())
            })
        })
    }

    /// Aligns a batch of paired-end reads.
    ///
    /// # Arguments
    /// - `records` - A slice of FASTQ read pairs.
    ///
    /// # Returns
    /// - An iterator over tuples of aligned SAM records for each mate.
    pub fn align_read_pairs<'a>(
        &'a self,
        records: &'a [(fastq::Record, fastq::Record)],
    ) -> impl ParallelIterator<Item = (Vec<sam::Record>, Vec<sam::Record>)> + 'a {
        let chunk_size = self.get_chunk_size(records.len());
        self.init_thread_pool().install(|| {
            records.par_chunks(chunk_size).flat_map_iter(|chunk| {
                let mut aligner = self.clone_inner_aligner();
                chunk
                    .iter()
                    .map(move |(fq1, fq2)| aligner.align_read_pair(fq1, fq2).unwrap())
            })
        })
    }

    /// This will return a copy of the aligner with the same reference index.
    /// Returning a new aligner is necessary for multithreaded alignment, as the
    /// underlying foreign aligner is not thread-safe.
    pub fn clone_inner_aligner(&self) -> InnerAligner {
        InnerAligner::new(self.index.as_ref())
    }

    fn init_thread_pool(&self) -> rayon::ThreadPool {
        rayon::ThreadPoolBuilder::new()
            .num_threads(self.opts.num_threads as usize)
            .build()
            .unwrap()
    }

    fn get_chunk_size(&self, n_records: usize) -> usize {
        (n_records + self.opts.num_threads as usize - 1) / self.opts.num_threads as usize
    }
}

fn parse_sam_to_records<'a>(
    aln_buf: &'a [u8],
    sam_buf: &'a mut Vec<u8>,
    name: &'a [u8],
) -> impl Iterator<Item = sam::Record> + 'a {
    aln_buf.split(|c| *c == b'\n').filter_map(|slc| {
        if !slc.is_empty() {
            sam_buf.clear();
            sam_buf.extend_from_slice(name);
            sam_buf.extend_from_slice(slc);
            Some(sam_buf.as_slice().try_into().unwrap())
        } else {
            None
        }
    })
}

/// Copy a fastq record into a buffer
fn copy_fq_to_buf(buf: &mut Vec<u8>, fq: &fastq::Record) -> Result<()> {
    buf.clear();
    buf.extend_from_slice(b"@a\n"); // make a fake name to avoid the copy of fq.name
    buf.extend_from_slice(fq.sequence());
    buf.extend_from_slice(b"\n+\n");
    buf.extend_from_slice(fq.quality_scores());
    buf.push(b'\n');
    buf.push(b'\0');
    Ok(())
}

/// Produces a header from the genome reference directory by looking up the contig names and
/// lengths and formatting them properly
fn generate_header<P: AsRef<Path>>(genome_path: P) -> Result<Header> {
    let mut header = Header::builder();

    let contig_names_path = genome_path.as_ref().join(Path::new("chrName.txt"));
    let contig_names = get_lines(&contig_names_path);

    let contig_lengths_path = genome_path.as_ref().join(Path::new("chrLength.txt"));
    let contig_lengths = get_lines(&contig_lengths_path);
    for (contig_name, len) in contig_names.iter().zip(contig_lengths.into_iter()) {
        header = header.add_reference_sequence(
            contig_name.as_bytes(),
            Map::<ReferenceSequence>::new(std::num::NonZeroUsize::try_from(
                len.parse::<usize>().unwrap(),
            )?),
        );
    }
    Ok(header.build())
}

fn get_lines(path: &Path) -> Vec<String> {
    let file = match File::open(path) {
        Err(error) => panic!("error: {}: {}", path.display(), error),
        Ok(file) => file,
    };
    BufReader::new(file).lines().map(Result::unwrap).collect()
}

#[cfg(test)]
mod test {
    use super::*;

    use fastq::record::Definition;
    use flate2::read::GzDecoder;

    fn make_fq(name: &[u8], seq: &[u8], qual: &[u8]) -> fastq::Record {
        fastq::Record::new(Definition::new(name, ""), seq, qual)
    }

    const ERCC_REF: &str = "test/ercc92-1.2.0/star/";

    const NAME: &[u8] = b"NAME";
    const ERCC_READ_1: &[u8] = b"GCATCCAGACCGTCGGCTGATCGTGGTTTTACTAGGCTAGACTAGCGTACGAGCACTATGGTCAGTAATTCCTGGAGGAATAGGTACCAAGAAAAAAACG";
    const ERCC_QUAL_1: &[u8] = b"????????????????????????????????????????????????????????????????????????????????????????????????????";

    const ERCC_READ_2: &[u8] = b"GGAGACGAATTGCCAGAATTATTAACTGCGCAGTTAGGGCAGCGTCTGAGGAAGTTTGCTGCGGTTTCGCCTTGACCGCGGGAAGGAGACATAACGATAG";
    const ERCC_QUAL_2: &[u8] = b"????????????????????????????????????????????????????????????????????????????????????????????????????";

    const ERCC_READ_3: &[u8] = b"AACTTAATGGACGGG";
    const ERCC_QUAL_3: &[u8] = b"???????????????";

    const ERCC_READ_4: &[u8] = b"AATCCACTCAATAAATCTAAAAAC";
    const ERCC_QUAL_4: &[u8] = b"????????????????????????";

    const ERCC_PARA300_PATH: &str = "test/test_ercc_reads.fastq.gz";

    #[test]
    fn test_header() {
        let opts = StarOpts::new(ERCC_REF);
        let aligner = StarAligner::new(opts).unwrap();

        let mut writer = sam::io::Writer::new(Vec::new());
        writer.write_header(&aligner.get_header()).unwrap();
        let header_string = String::from_utf8(writer.get_ref().clone()).unwrap();
        assert!(header_string.starts_with("@SQ\tSN:ERCC-00002\tLN:1061\n"));
        assert!(header_string.ends_with("@SQ\tSN:ERCC-00171\tLN:505\n"));
    }

    #[test]
    fn test_empty_tiny_reads() {
        let opts = StarOpts::new(ERCC_REF);
        let mut aligner = StarAligner::new(opts).unwrap().clone_inner_aligner();

        let recs = aligner.align_read(&make_fq(b"a", b"b", b"?")).unwrap();
        println!("{:?}", recs);

        let (recs1, recs2) = aligner.align_read_pair(
            &make_fq(b"a", b"A", b"?"),
            &make_fq(b"a", b"C", b"?"),
        )
        .unwrap();
        println!("{:?}, {:?}", recs1, recs2);
    }

    #[test]
    fn test_ercc_align() {
        let opts = StarOpts::new(ERCC_REF);
        let aligner = StarAligner::new(opts).unwrap();
        let header = aligner.get_header().clone();

        let fq1 = make_fq(b"ERCC1", ERCC_READ_1, ERCC_QUAL_1);
        let fq2 = make_fq(b"ERCC2", ERCC_READ_2, ERCC_QUAL_2);
        let fq3 = make_fq(b"ERCC3", ERCC_READ_3, ERCC_QUAL_3);
        let fq4 = make_fq(b"ERCC4", ERCC_READ_4, ERCC_QUAL_4);

        let recs = aligner.clone_inner_aligner().align_read(&fq1).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].name().unwrap(), "ERCC1");
        assert_eq!(recs[0].alignment_start().unwrap().unwrap().get(), 51);
        assert_eq!(recs[0].reference_sequence_id(&header).unwrap().unwrap(), 0);
        println!("{:?}", recs);

        let recs = aligner.clone_inner_aligner().align_read(&fq2).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].reference_sequence_id(&header).unwrap().unwrap(), 0);
        assert_eq!(recs[0].alignment_start().unwrap().unwrap().get(), 501);
        println!("{:?}", recs);

        let recs = aligner.clone_inner_aligner().align_read(&fq3).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].flags().unwrap().bits(), 0);
        assert_eq!(recs[0].reference_sequence_id(&header).unwrap().unwrap(), 39);
        assert_eq!(recs[0].alignment_start().unwrap().unwrap().get(), 28);
        assert_eq!(recs[0].mapping_quality().unwrap().unwrap().get(), 3);
        assert_eq!(recs[1].flags().unwrap().bits(), 0x110); // REVERSE,SECONDARY
        assert_eq!(recs[1].reference_sequence_id(&header).unwrap().unwrap(), 72);
        assert_eq!(recs[1].alignment_start().unwrap().unwrap().get(), 554);
        assert_eq!(recs[1].mapping_quality().unwrap().unwrap().get(), 3);
        println!("{:?}", recs);

        let recs = aligner.clone_inner_aligner().align_read(&fq4).unwrap();
        println!("{:?}", recs);
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].flags().unwrap().bits(), 0);
        assert_eq!(recs[0].reference_sequence_id(&header).unwrap().unwrap(), 72);
        assert_eq!(recs[0].alignment_start().unwrap().unwrap().get(), 493);
        assert_eq!(recs[0].mapping_quality().unwrap().unwrap().get(), 3);
        assert_eq!(recs[1].flags().unwrap().bits(), 0x100); // SECONDARY
        assert_eq!(recs[1].reference_sequence_id(&header).unwrap().unwrap(), 72);
        assert_eq!(recs[1].alignment_start().unwrap().unwrap().get(), 608);
        assert_eq!(recs[1].mapping_quality().unwrap().unwrap().get(), 3);
    }

    #[test]
    fn test_transcriptome_min_score() {
        let mut opts = StarOpts::new(ERCC_REF);
        opts.out_filter_score_min = 20;
        let aligner = StarAligner::new(opts).unwrap();
        let header = aligner.get_header().clone();

        let fq = make_fq(NAME, ERCC_READ_3, ERCC_QUAL_3);

        let recs = aligner.clone_inner_aligner().align_read(&fq).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].flags().unwrap().bits(), 4); // UNMAP
        assert!(recs[0].reference_sequence_id(&header).is_none());
        assert!(recs[0].alignment_start().is_none());
        assert_eq!(recs[0].mapping_quality().unwrap().unwrap().get(), 0);
        println!("{:?}", recs);
    }

    #[test]
    fn test_multithreads_align() {
        let mut fq_reader = File::open(ERCC_PARA300_PATH)
            .map(GzDecoder::new)
            .map(BufReader::new)
            .map(fastq::io::Reader::new)
            .unwrap();
        let records: Vec<_> = fq_reader.records().map(Result::unwrap).collect();

        // Single end
        let opts = StarOpts::new(ERCC_REF);
        let aligner = StarAligner::new(opts).unwrap().with_num_threads(1);
        let res1 = aligner.align_reads(&records).collect::<Vec<_>>();

        let aligner = aligner.with_num_threads(8);
        let res2 = aligner.align_reads(&records).collect::<Vec<_>>();
        assert_eq!(res1, res2);

        // Paired end
        let records = records
            .chunks(2)
            .map(|chunk| {
                let (r1, r2) = chunk.split_at(1);
                (r1[0].clone(), r2[0].clone())
            })
            .collect::<Vec<_>>();

        let aligner = aligner.with_num_threads(1);
        let res1 = aligner.align_read_pairs(&records).collect::<Vec<_>>();

        let aligner = aligner.with_num_threads(8);
        let res2 = aligner.align_read_pairs(&records).collect::<Vec<_>>();
        assert_eq!(res1, res2);
    }
}
