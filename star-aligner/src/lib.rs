use anyhow::{bail, Result};
use noodles::fastq::Record as FastqRecord;
use noodles::sam::{
    header::record::value::{map::ReferenceSequence, Map},
    Header, Record as SamRecord,
};
use noodles::{fastq, sam};
use rayon::prelude::*;
use std::{
    ffi::{c_char, c_int, CStr, CString},
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    sync::Arc,
    vec,
};

/// StarOpts contains the parameters which will be used for the STAR aligner.
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
    pub fn new<P: AsRef<Path>>(genome_dir: P) -> Self {
        Self {
            genome_dir: genome_dir.as_ref().to_path_buf(),
            num_threads: 8,
            sjdb_overhang: 100,
            out_filter_score_min: 0,
        }
    }

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

/// A wrapper around the STAR aligner. This is solely used to ensure that the
/// aligner is properly deallocated when it goes out of scope.
struct _AlignerWrapper(*mut star_sys::Aligner);

impl _AlignerWrapper {
    fn new(index: &StarIndex) -> Self {
        unsafe { Self(star_sys::init_aligner_from_ref(index.index)) }
    }
}

impl Drop for _AlignerWrapper {
    fn drop(&mut self) {
        unsafe { star_sys::destroy_aligner(self.0) };
    }
}

pub struct StarAligner {
    index: Arc<StarIndex>,
    opts: StarOpts,
}

/// StarAligner aligns single reads or read-pairs to the reference it is initialized with, and returns
/// rust_htslib Record objects
impl StarAligner {
    /// Load the reference index into memory based on the given settings.
    /// Should only be called once for a given reference. Create aligners
    /// that use this reference by calling `get_aligner`.  The reference
    /// index will be free'd when all `Aligner`s that use this `StarReference`
    /// have been dropped.
    pub fn new(opts: StarOpts) -> Result<Self> {
        let header = generate_header(opts.genome_dir.as_path())?;

        let c_args = opts.to_c_args();
        let length = c_args.len() as c_int;

        let index = unsafe { star_sys::init_star_ref(length, c_args.as_ptr()) };
        let index = StarIndex { index, header };

        Ok(Self {
            index: Arc::new(index),
            opts,
        })
    }

    /// Set the number of threads to use for alignment. Default is 8.
    pub fn with_num_threads(mut self, n: u16) -> Self {
        self.opts.num_threads = n;
        self
    }

    pub fn get_header(&self) -> &Header {
        &self.index.header
    }

    /// Align a chunk of reads
    pub fn align_reads<'a>(
        &'a self,
        records: &'a [fastq::Record],
    ) -> impl ParallelIterator<Item = Vec<SamRecord>> + 'a {
        let chunk_size = self.get_chunk_size(records.len());
        self.init_thread_pool().install(|| {
            records.par_chunks(chunk_size).flat_map(|chunk| {
                let mut aligner = self.get_aligner();
                let mut fq_buf = Vec::new();
                let mut sam_buf = Vec::new();
                chunk
                    .iter()
                    .map(|fq| align_read(&mut aligner, fq, &mut fq_buf, &mut sam_buf).unwrap())
                    .collect::<Vec<_>>()
            })
        })
    }

    /// Align a chunk of read pairs
    pub fn align_read_pairs<'a>(
        &'a self,
        records: &'a [(fastq::Record, fastq::Record)],
    ) -> impl ParallelIterator<Item = (Vec<SamRecord>, Vec<SamRecord>)> + 'a {
        let chunk_size = self.get_chunk_size(records.len());
        self.init_thread_pool().install(|| {
            records.par_chunks(chunk_size).flat_map(|chunk| {
                let mut aligner = self.get_aligner();
                let mut fq_buf1 = Vec::new();
                let mut fq_buf2 = Vec::new();
                let mut sam_buf = Vec::new();
                chunk
                    .iter()
                    .map(|(fq1, fq2)| {
                        align_read_pair(
                            &mut aligner,
                            fq1,
                            fq2,
                            &mut fq_buf1,
                            &mut fq_buf2,
                            &mut sam_buf,
                        )
                        .unwrap()
                    })
                    .collect::<Vec<_>>()
            })
        })
    }

    /// This will return a copy of the aligner with the same reference index.
    /// Returning a new aligner is necessary for multithreaded alignment, as the
    /// underlying foreign aligner is not thread-safe.
    fn get_aligner(&self) -> _AlignerWrapper {
        _AlignerWrapper::new(self.index.as_ref())
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

/// Align a single read. Note that the aligner is mutably borrowed.
/// We did this on purpose to ensure this function cannot be called concurrently.
fn align_read(
    aligner: &mut _AlignerWrapper,
    fq: &FastqRecord,
    fq_buf: &mut Vec<u8>,
    sam_buf: &mut Vec<u8>,
) -> Result<Vec<SamRecord>> {
    // Copy the fastq record into a buffer
    copy_fq_to_buf(fq_buf, fq)?;

    let fq_cstr = CStr::from_bytes_with_nul(fq_buf)?;
    let res: *const c_char = unsafe { star_sys::align_read(aligner.0, fq_cstr.as_ptr()) };
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
fn align_read_pair(
    aligner: &mut _AlignerWrapper,
    fq1: &FastqRecord,
    fq2: &FastqRecord,
    fq_buf1: &mut Vec<u8>,
    fq_buf2: &mut Vec<u8>,
    sam_buf: &mut Vec<u8>,
) -> Result<(Vec<SamRecord>, Vec<SamRecord>)> {
    // Copy the fastq records into buffers
    copy_fq_to_buf(fq_buf1, fq1)?;
    copy_fq_to_buf(fq_buf2, fq2)?;
    let fq1_str = CStr::from_bytes_with_nul(fq_buf1)?;
    let fq2_str = CStr::from_bytes_with_nul(fq_buf2)?;

    let res: *const c_char =
        unsafe { star_sys::align_read_pair(aligner.0, fq1_str.as_ptr(), fq2_str.as_ptr()) };
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
fn copy_fq_to_buf(buf: &mut Vec<u8>, fq: &FastqRecord) -> Result<()> {
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

    fn make_fq(name: &[u8], seq: &[u8], qual: &[u8]) -> FastqRecord {
        FastqRecord::new(Definition::new(name, ""), seq, qual)
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
        writer.write_header(aligner.get_header()).unwrap();
        let header_string = String::from_utf8(writer.get_ref().clone()).unwrap();
        assert!(header_string.starts_with("@SQ\tSN:ERCC-00002\tLN:1061\n"));
        assert!(header_string.ends_with("@SQ\tSN:ERCC-00171\tLN:505\n"));
    }

    #[test]
    fn test_empty_tiny_reads() {
        let opts = StarOpts::new(ERCC_REF);
        let mut aligner = StarAligner::new(opts).unwrap().get_aligner();

        let mut fq_buf1 = Vec::new();
        let mut fq_buf2 = Vec::new();
        let mut sam_buf = Vec::new();
        let recs = align_read(
            &mut aligner,
            &make_fq(b"a", b"b", b"?"),
            &mut fq_buf1,
            &mut sam_buf,
        )
        .unwrap();
        println!("{:?}", recs);

        let (recs1, recs2) = align_read_pair(
            &mut aligner,
            &make_fq(b"a", b"A", b"?"),
            &make_fq(b"a", b"C", b"?"),
            &mut fq_buf1,
            &mut fq_buf2,
            &mut sam_buf,
        )
        .unwrap();
        println!("{:?}, {:?}", recs1, recs2);
    }

    #[test]
    fn test_ercc_align() {
        let opts = StarOpts::new(ERCC_REF);
        let aligner = StarAligner::new(opts).unwrap();
        let header = aligner.get_header().clone();

        let fq1 = make_fq(NAME, ERCC_READ_1, ERCC_QUAL_1);
        let fq2 = make_fq(NAME, ERCC_READ_2, ERCC_QUAL_2);
        let fq3 = make_fq(NAME, ERCC_READ_3, ERCC_QUAL_3);
        let fq4 = make_fq(NAME, ERCC_READ_4, ERCC_QUAL_4);

        let mut fq_buf = Vec::new();
        let mut sam_buf = Vec::new();

        let recs = align_read(&mut aligner.get_aligner(), &fq1, &mut fq_buf, &mut sam_buf).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].alignment_start().unwrap().unwrap().get(), 51);
        assert_eq!(recs[0].reference_sequence_id(&header).unwrap().unwrap(), 0);
        println!("{:?}", recs);

        let recs = align_read(&mut aligner.get_aligner(), &fq2, &mut fq_buf, &mut sam_buf).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].reference_sequence_id(&header).unwrap().unwrap(), 0);
        assert_eq!(recs[0].alignment_start().unwrap().unwrap().get(), 501);
        println!("{:?}", recs);

        let recs = align_read(&mut aligner.get_aligner(), &fq3, &mut fq_buf, &mut sam_buf).unwrap();
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

        let recs = align_read(&mut aligner.get_aligner(), &fq4, &mut fq_buf, &mut sam_buf).unwrap();
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
        let mut fq_buf = Vec::new();
        let mut sam_buf = Vec::new();

        let recs = align_read(&mut aligner.get_aligner(), &fq, &mut fq_buf, &mut sam_buf).unwrap();
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
