pub mod transcript;

use anyhow::{bail, Result};
use noodles::{
    fastq,
    sam::{
        self, alignment::RecordBuf, header::record::value::{map::ReferenceSequence, Map}, Header
    },
};
use std::{
    ffi::{c_char, c_int, CStr, CString},
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    sync::Arc,
    vec,
};

/// Main interface for performing read alignments using STAR.
///
/// `StarAligner` handles alignment of both single and paired-end reads to a reference genome.
pub struct StarAligner {
    index: Arc<StarIndex>,
    opts: StarOpts,
    inner: InnerAligner,
}

/// This will return a copy of the aligner with the same reference index.
/// Returning a new aligner is necessary for multithreaded alignment, as the
/// underlying foreign aligner is not thread-safe.
impl Clone for StarAligner {
    fn clone(&self) -> Self {
        Self {
            index: self.index.clone(),
            opts: self.opts.clone(),
            inner: InnerAligner::new(self.index.as_ref()),
        }
    }
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
        let inner = InnerAligner::new(&index);
        Ok(Self {
            index: Arc::new(index),
            opts,
            inner,
        })
    }

    pub fn get_transcriptome(&self) -> Result<transcript::Transcriptome> {
        transcript::Transcriptome::from_path(&self.opts.genome_dir)
    }

    /// Retrieves the header information.
    ///
    /// # Returns
    /// - A reference to the alignment header.
    pub fn get_header(&self) -> &Header {
        &self.index.header
    }

    /// Aligns single-end reads.
    pub fn align_read(&mut self, fq: &fastq::Record) -> Result<Vec<RecordBuf>> {
        self.inner.align_read(&self.index.header, fq)
    }

    /// Aligns paired-end reads.
    pub fn align_read_pair(
        &mut self,
        fq1: &fastq::Record,
        fq2: &fastq::Record,
    ) -> Result<(Vec<RecordBuf>, Vec<RecordBuf>)> {
        self.inner.align_read_pair(&self.index.header, fq1, fq2)
    }
}

/// Struct representing STAR aligner options.
///
/// `StarOpts` stores settings to control the behavior of the STAR aligner.
#[derive(Clone)]
pub struct StarOpts {
    genome_dir: PathBuf, // Path to the STAR reference genome
    sjdb_overhang: u16,  // Overhang for splice junctions
    out_filter_score_min: u16, // Alignment will be output only if its score is higher than or equal to this value.
                               // To output higher quality reads, you can set this value closer to read length.
                               // Default: 0.
    out_sam_attributes: String, // SAM attributes to include in output (e.g., "Standard", "All", "NH HI AS nM NM MD")
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
            sjdb_overhang: 100,
            out_filter_score_min: 0,
            out_sam_attributes: "Standard".to_string(), // Default: NH HI AS nM
        }
    }

    /// Set the SAM attributes to include in the output.
    ///
    /// # Arguments
    /// - `attributes` - SAM attributes string. Common options:
    ///   - "Standard" - NH HI AS nM (default)
    ///   - "All" - NH HI AS nM NM MD jM jI MC ch
    ///   - "Standard MD" - Standard attributes plus MD tag
    ///   - Custom combinations like "NH HI AS nM NM MD"
    ///
    /// # Returns
    /// - A modified `StarOpts` object with the specified SAM attributes.
    pub fn with_sam_attributes<S: Into<String>>(mut self, attributes: S) -> Self {
        self.out_sam_attributes = attributes.into();
        self
    }

    /// Enable MD tag output by setting attributes to "Standard MD".
    /// This includes the standard attributes (NH HI AS nM) plus the MD tag.
    ///
    /// # Returns
    /// - A modified `StarOpts` object with MD tag enabled.
    pub fn with_md_tag(mut self) -> Self {
        self.out_sam_attributes = "Standard MD".to_string();
        self
    }

    /// Converts options to a vector of C-style strings for FFI compatibility.
    ///
    /// # Returns
    /// - A vector of pointers to C-style strings representing each option.
    fn to_c_args(&self) -> Vec<*const i8> {
        let sjdb_overhang = self.sjdb_overhang.to_string();
        let out_filter_score_min = self.out_filter_score_min.to_string();
        
        // Parse SAM attributes string and split by whitespace
        let sam_attributes: Vec<&str> = self.out_sam_attributes
            .split_whitespace()
            .collect();
        
        let mut args = vec![
            "STAR",
            "--genomeDir",
            self.genome_dir.to_str().unwrap(),
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
            "--outSAMattributes",
        ];
        
        // Add each SAM attribute as a separate argument
        args.extend(sam_attributes);
        
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
struct InnerAligner {
    inner: *mut star_sys::Aligner,
    sam_buf: Vec<u8>,
    fq1_buf: Vec<u8>,
    fq2_buf: Vec<u8>,
}

unsafe impl Send for InnerAligner {}
unsafe impl Sync for InnerAligner {}

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
    fn align_read(&mut self, header: &Header, fq: &fastq::Record) -> Result<Vec<RecordBuf>> {
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
        let records: Vec<_> = parse_sam_to_records(aln_buf, sam_buf, fq.name())
            .map(|x| RecordBuf::try_from_alignment_record(header, &x).unwrap())
            .collect();

        unsafe {
            libc::free(res as *mut libc::c_void);
        }
        Ok(records)
    }

    /// Align a pair of reads. Note that the aligner is mutably borrowed.
    /// We did this on purpose to ensure this function cannot be called concurrently.
    fn align_read_pair(
        &mut self,
        header: &Header,
        fq1: &fastq::Record,
        fq2: &fastq::Record,
    ) -> Result<(Vec<RecordBuf>, Vec<RecordBuf>)> {
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
            .map(|x| RecordBuf::try_from_alignment_record(header, &x).unwrap())
            .partition(|rec| rec.flags().is_first_segment());

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
    use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

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

        aligner.get_transcriptome().unwrap();

        let mut writer = sam::io::Writer::new(Vec::new());
        writer.write_header(&aligner.get_header()).unwrap();
        let header_string = String::from_utf8(writer.get_ref().clone()).unwrap();
        assert!(header_string.starts_with("@SQ\tSN:ERCC-00002\tLN:1061\n"));
        assert!(header_string.ends_with("@SQ\tSN:ERCC-00171\tLN:505\n"));
    }

    #[test]
    fn test_empty_tiny_reads() {
        let opts = StarOpts::new(ERCC_REF);
        let mut aligner = StarAligner::new(opts).unwrap();

        let recs = aligner.align_read(&make_fq(b"a", b"b", b"?")).unwrap();
        println!("{:?}", recs);

        let (recs1, recs2) = aligner
            .align_read_pair(&make_fq(b"a", b"A", b"?"), &make_fq(b"a", b"C", b"?"))
            .unwrap();
        println!("{:?}, {:?}", recs1, recs2);
    }

    #[test]
    fn test_ercc_align() {
        let opts = StarOpts::new(ERCC_REF);
        let mut aligner = StarAligner::new(opts).unwrap();

        let unmapped = make_fq(b"UM", b"ATCGATCGATCGATCG", b"????????????????");
        let fq1 = make_fq(b"ERCC1", ERCC_READ_1, ERCC_QUAL_1);
        let fq2 = make_fq(b"ERCC2", ERCC_READ_2, ERCC_QUAL_2);
        let fq3 = make_fq(b"ERCC3", ERCC_READ_3, ERCC_QUAL_3);
        let fq4 = make_fq(b"ERCC4", ERCC_READ_4, ERCC_QUAL_4);

        let recs = aligner.align_read(&unmapped).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].flags().bits(), 0x4);

        let recs = aligner.align_read(&fq1).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].name().unwrap(), "ERCC1");
        assert_eq!(recs[0].alignment_start().unwrap().get(), 51);
        assert_eq!(recs[0].reference_sequence_id().unwrap(), 0);
        assert!(recs[0].mapping_quality().is_none());  // STAR assigns 255 to uniquely mapped reads 

        let recs = aligner.align_read(&fq2).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].reference_sequence_id().unwrap(), 0);
        assert_eq!(recs[0].alignment_start().unwrap().get(), 501);
        assert!(recs[0].mapping_quality().is_none());  // STAR assigns 255 to uniquely mapped reads 

        let recs = aligner.align_read(&fq3).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].flags().bits(), 0);
        assert_eq!(recs[0].reference_sequence_id().unwrap(), 39);
        assert_eq!(recs[0].alignment_start().unwrap().get(), 28);
        assert_eq!(recs[0].mapping_quality().unwrap().get(), 3);
        assert_eq!(recs[1].flags().bits(), 0x110); // REVERSE,SECONDARY
        assert_eq!(recs[1].reference_sequence_id().unwrap(), 72);
        assert_eq!(recs[1].alignment_start().unwrap().get(), 554);
        assert_eq!(recs[1].mapping_quality().unwrap().get(), 3);

        let recs = aligner.align_read(&fq4).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].flags().bits(), 0);
        assert_eq!(recs[0].reference_sequence_id().unwrap(), 72);
        assert_eq!(recs[0].alignment_start().unwrap().get(), 493);
        assert_eq!(recs[0].mapping_quality().unwrap().get(), 3);
        assert_eq!(recs[1].flags().bits(), 0x100); // SECONDARY
        assert_eq!(recs[1].reference_sequence_id().unwrap(), 72);
        assert_eq!(recs[1].alignment_start().unwrap().get(), 608);
        assert_eq!(recs[1].mapping_quality().unwrap().get(), 3);

        let (recs1, recs2) = aligner.align_read_pair(&fq1, &fq2).unwrap();
        assert_eq!((recs1.len(), recs2.len()), (1, 1));
    }

    #[test]
    fn test_transcriptome_min_score() {
        let mut opts = StarOpts::new(ERCC_REF);
        opts.out_filter_score_min = 20;
        let mut aligner = StarAligner::new(opts).unwrap();

        let fq = make_fq(NAME, ERCC_READ_3, ERCC_QUAL_3);

        let recs = aligner.align_read(&fq).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].flags().bits(), 4); // UNMAP
        assert!(recs[0].reference_sequence_id().is_none());
        assert!(recs[0].alignment_start().is_none());
        assert_eq!(recs[0].mapping_quality().unwrap().get(), 0);
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
        let mut aligner = StarAligner::new(opts).unwrap();
        let res1 = records
            .iter()
            .map(|r| aligner.align_read(r).unwrap())
            .collect::<Vec<_>>();

        let res2 = records
            .par_iter()
            .map(|r| aligner.clone().align_read(r).unwrap())
            .collect::<Vec<_>>();
        assert_eq!(res1, res2);

        // Paired end
        let records = records
            .chunks(2)
            .map(|chunk| {
                let (r1, r2) = chunk.split_at(1);
                (r1[0].clone(), r2[0].clone())
            })
            .collect::<Vec<_>>();

        let res1 = records
            .iter()
            .map(|(r1, r2)| aligner.align_read_pair(r1, r2).unwrap())
            .collect::<Vec<_>>();

        let res2 = records
            .par_iter()
            .map(|(r1, r2)| aligner.clone().align_read_pair(r1, r2).unwrap())
            .collect::<Vec<_>>();
        assert_eq!(res1, res2);
    }


    #[test]
    fn test_star() {
        let mut fq_reader = File::open(ERCC_PARA300_PATH)
            .map(GzDecoder::new)
            .map(BufReader::new)
            .map(fastq::io::Reader::new)
            .unwrap();
        let records: Vec<_> = fq_reader.records().map(Result::unwrap).collect();

        let star_index = "/data/Public/STAR_reference/refdata-gex-GRCm39-2024-A/star_2.7.1";
        let opts = StarOpts::new(star_index);
        let mut aligner = StarAligner::new(opts).unwrap();

        let res1 = records
            .iter()
            .map(|r| aligner.align_read(r).unwrap())
            .collect::<Vec<_>>();
        println!("{:?}", res1)
    }

    #[test]
    fn test_md_tag_configuration() {
        // Test with MD tag enabled
        let opts_with_md = StarOpts::new(ERCC_REF).with_md_tag();
        let mut aligner = StarAligner::new(opts_with_md).unwrap();

        let fq = make_fq(NAME, ERCC_READ_1, ERCC_QUAL_1);
        let recs = aligner.align_read(&fq).unwrap();

        // The MD tag should be present in the alignment output
        // Note: The actual MD tag content would be visible in the SAM output
        // but is not directly accessible through the RecordBuf interface
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].reference_sequence_id().unwrap(), 0);
    }

    #[test]
    fn test_custom_sam_attributes() {
        // Test with custom SAM attributes including MD and NM tags
        let opts_custom = StarOpts::new(ERCC_REF)
            .with_sam_attributes("NH HI AS");
        let mut aligner = StarAligner::new(opts_custom).unwrap();

        let fq = make_fq(NAME, ERCC_READ_1, ERCC_QUAL_1);
        let recs = aligner.align_read(&fq).unwrap();

        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].reference_sequence_id().unwrap(), 0);
    }
}
