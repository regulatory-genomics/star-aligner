use std::{collections::HashMap, ffi::{c_char, c_int, CStr, CString}, fs::File, io::{BufRead, BufReader}, path::{Path, PathBuf}, sync::Arc, vec};
use noodles::sam::{header::record::value::{Map, map::ReferenceSequence}, Header};
use anyhow::{bail, Result};
use noodles::{fastq, sam};

struct InnerStarIndex {
    index: *const star_sys::StarRef,
    opts: StarOpts,
    header: Header,
}

impl Drop for InnerStarIndex {
    fn drop(&mut self) {
        unsafe {
            star_sys::destroy_ref(self.index);
        }
    }
}

#[derive(Clone)]
pub struct StarIndex(Arc<InnerStarIndex>);

impl StarIndex {
    /// Load the reference index into memory based on the given settings.
    /// Should only be called once for a given reference. Create aligners
    /// that use this reference by calling `get_aligner`.  The reference
    /// index will be free'd when all `Aligner`s that use this `StarReference`
    /// have been dropped.
    pub fn load(opts: StarOpts) -> Result<StarIndex> {
        let header = generate_header(opts.genome_dir.as_path())?;

        let c_args = opts.to_c_args();
        let length = c_args.len() as c_int;
        let index = unsafe { star_sys::init_star_ref(length, c_args.as_ptr()) };

        let inner = InnerStarIndex {
            index,
            opts,
            header,
        };

        Ok(StarIndex(Arc::new(inner)))
    }

    pub fn header(&self) -> &Header {
        &self.0.header
    }
}


/// StarOpts contains the parameters which will be used for the STAR aligner.
/// Currently the array of argument strings is passed directly
#[derive(Clone)]
pub struct StarOpts {
    genome_dir: PathBuf,  // Path to the STAR reference genome
    num_threads: u16,     // Number of threads to use
    sjdb_overhang: u16,   // Overhang for splice junctions
    out_filter_score_min: u16,  // Alignment will be output only if its score is higher than or equal to this value.
                               // To output higher quality reads, you can set this value closer to read length.
                               // Default: 0.
}

impl StarOpts {
    fn new<P: AsRef<Path>>(genome_dir: P) -> Self {
        Self {
            genome_dir: genome_dir.as_ref().to_path_buf(),
            num_threads: 1,
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
            "--genomeDir", self.genome_dir.to_str().unwrap(),
            "--runThreadN", num_threads.as_str(),
            "--sjdbOverhang", sjdb_overhang.as_str(),
            "--outFilterScoreMin", out_filter_score_min.as_str(),
            "--readNameSeparator", "space",
            "--outSAMunmapped", "Within", "KeepPairs",
            "--outSAMtype", "SAM",
            "--outStd", "SAM",
            "--outSAMorder", "PairedKeepInputOrder",
        ];
        args.into_iter().map(|s| CString::new(s).unwrap().into_raw() as *const i8).collect()
    }

}

/// StarAligner aligns single reads or read-pairs to the reference it is initialized with, and returns
/// rust_htslib Record objects
pub struct StarAligner {
    aligner: *mut star_sys::Aligner,
    index: StarIndex,
    sam_buf: Vec<u8>,
    aln_buf: Vec<u8>,
    fastq1: Vec<u8>,
    fastq2: Vec<u8>,
}

impl Drop for StarAligner {
    fn drop(&mut self) {
        unsafe { star_sys::destroy_aligner(self.aligner) };
    }
}

unsafe impl Send for StarAligner {}

impl StarAligner {
    pub fn new(index: StarIndex) -> Self {
        let aligner = unsafe { star_sys::init_aligner_from_ref(index.0.as_ref().index) };
        Self {
            aligner,
            index,
            sam_buf: Vec::new(),
            aln_buf: Vec::new(),
            fastq1: Vec::new(),
            fastq2: Vec::new(),
        }
    }

    pub fn get_header(&self) -> &Header {
        &self.index.0.header
    }

    fn align_read(&mut self, name: &[u8], seq: &[u8], qual: &[u8]) -> Vec<sam::Record> {
        Self::prepare_fastq(&mut self.fastq1, name, seq, qual);

        align_read_rust(self.aligner, self.fastq1.as_slice(), &mut self.aln_buf).unwrap();
        self.parse_sam_to_records(name)
    }

    pub fn align_read_pair(
        &mut self,
        name: &[u8],
        read1: &[u8],
        qual1: &[u8],
        read2: &[u8],
        qual2: &[u8],
    ) -> (Vec<sam::Record>, Vec<sam::Record>) {
        Self::prepare_fastq(&mut self.fastq1, name, read1, qual1);
        Self::prepare_fastq(&mut self.fastq2, name, read2, qual2);
        align_read_pair_rust(
            self.aligner,
            self.fastq1.as_slice(),
            self.fastq2.as_slice(),
            &mut self.aln_buf,
        )
        .unwrap();
        let full_vec = self.parse_sam_to_records(name);

        // Partition the records into first mate and second mate
        let mut first_vec: Vec<sam::Record> = Vec::new();
        let mut second_vec: Vec<sam::Record> = Vec::new();
        for rec in full_vec {
            if rec.flags().unwrap().is_first_segment() {
                first_vec.push(rec);
            } else {
                second_vec.push(rec);
            }
        }
        (first_vec, second_vec)
    }


    /// Prepare a read and qual for processing by formatting as FASTQ
    fn prepare_fastq(buf: &mut Vec<u8>, _name: &[u8], read: &[u8], qual: &[u8]) {
        buf.clear();
        // for now, we fixup the read name on the backend, avoid the copy for now
        buf.extend_from_slice(b"@a\n");
        buf.extend_from_slice(read);
        buf.extend_from_slice(b"\n+\n");
        buf.extend_from_slice(qual);
        buf.push(b'\n');
        buf.push(b'\0');
    }

    fn parse_sam_to_records(&mut self, name: &[u8]) -> Vec<sam::Record> {
        let mut records = Vec::new();
        for slc in self.aln_buf.split(|c| *c == b'\n') {
            if !slc.is_empty() {
                self.sam_buf.clear();
                self.sam_buf.extend_from_slice(name);
                self.sam_buf.extend_from_slice(slc);
                let record = self.sam_buf.as_slice().try_into().unwrap();
                records.push(record);
            }
        }

        records
    }
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
            Map::<ReferenceSequence>::new(std::num::NonZeroUsize::try_from(len.parse::<usize>().unwrap())?),
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

fn align_read_rust(al: *mut star_sys::Aligner, fastq: &[u8], aln_buf: &mut Vec<u8>) -> Result<()> {
    let fastq = CStr::from_bytes_with_nul(fastq)?;
    let res: *const c_char = unsafe { star_sys::align_read(al, fastq.as_ptr()) };
    if res.is_null() {
        bail!("STAR returned null alignment");
    }

    let cstr = unsafe { CStr::from_ptr(res) };
    aln_buf.clear();
    aln_buf.extend_from_slice(cstr.to_bytes());

    unsafe {
        libc::free(res as *mut libc::c_void);
    }
    Ok(())
}

fn align_read_pair_rust(
    al: *mut star_sys::Aligner,
    fastq1: &[u8],
    fastq2: &[u8],
    aln_buf: &mut Vec<u8>,
) -> Result<()> {
    let fastq1 = CStr::from_bytes_with_nul(fastq1)?;
    let fastq2 = CStr::from_bytes_with_nul(fastq2)?;
    let res: *const c_char =
        unsafe { star_sys::align_read_pair(al, fastq1.as_ptr(), fastq2.as_ptr()) };
    if res.is_null() {
        bail!("STAR returned null alignment");
    }

    let cstr = unsafe { CStr::from_ptr(res) };
    aln_buf.clear();
    aln_buf.extend_from_slice(cstr.to_bytes());

    unsafe {
        libc::free(res as *mut libc::c_void);
    }
    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;

    /// References to some commonly used reference genomes for testing purposes
    pub const DEFAULT_REF_1: &str = "/mnt/opt/refdata_cellranger/mm10-3.0.0/star";
    pub const DEFAULT_REF_2: &str = "/mnt/opt/refdata_cellranger/GRCh38-3.0.0/star";
    pub const DEFAULT_REF_3: &str = "/mnt/opt/refdata_cellranger/GRCh38-1.2.0/star";
    pub const DEFAULT_REF_4: &str = "//mnt/opt/refdata_cellranger/hg19_and_mm10-3.0.0/star";

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

    fn have_refs() -> bool {
        Path::new("/mnt/opt/refdata_cellranger").exists()
    }

     #[test]
    fn test_empty_tiny_reads() {
        let opts = StarOpts::new(ERCC_REF);
        let index = StarIndex::load(opts).unwrap();
        let mut aligner = StarAligner::new(index);

        let recs = aligner.align_read(b"b", b"A", b"?");
        println!("{:?}", recs);
        let (recs1, recs2) = aligner.align_read_pair(b"d", b"A", b"?", b"C", b"?");
        println!("{:?}, {:?}", recs1, recs2);
    }

    #[test]
    fn test_ercc_align() {
        let opts = StarOpts::new(ERCC_REF);
        let index = StarIndex::load(opts).unwrap();
        let mut aligner = StarAligner::new(index);
        let header = aligner.index.0.header.clone();

        let recs = aligner.align_read(NAME, ERCC_READ_1, ERCC_QUAL_1);
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].alignment_start().unwrap().unwrap().get(), 51);
        assert_eq!(recs[0].reference_sequence_id(&header).unwrap().unwrap(), 0);
        println!("{:?}", recs);

        let recs = aligner.align_read(NAME, ERCC_READ_2, ERCC_QUAL_2);
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].reference_sequence_id(&header).unwrap().unwrap(), 0);
        assert_eq!(recs[0].alignment_start().unwrap().unwrap().get(), 501);
        println!("{:?}", recs);

        let recs = aligner.align_read(NAME, ERCC_READ_3, ERCC_QUAL_3);
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

        let recs = aligner.align_read(NAME, ERCC_READ_4, ERCC_QUAL_4);
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
        let index = StarIndex::load(opts).unwrap();
        let mut aligner = StarAligner::new(index);
        let header = aligner.index.0.header.clone();

        let recs = aligner.align_read(NAME, ERCC_READ_3, ERCC_QUAL_3);
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].flags().unwrap().bits(), 4); // UNMAP
        assert!(recs[0].reference_sequence_id(&header).is_none());
        assert!(recs[0].alignment_start().is_none());
        assert_eq!(recs[0].mapping_quality().unwrap().unwrap().get(), 0);
        println!("{:?}", recs);
    }

    #[test]
    fn test_multithreaded_alignment() {
        let opts = StarOpts::new(ERCC_REF);
        let index = StarIndex::load(opts).unwrap();
        let mut aligner1 = StarAligner::new(index.clone());
        let mut aligner2 = StarAligner::new(index);

        let t1 = std::thread::spawn(move || {
            let header = aligner1.index.0.header.clone();
            for _ in 0..100000 {
                let recs = aligner1.align_read(NAME, ERCC_READ_1, ERCC_QUAL_1);
                assert_eq!(recs.len(), 1);
                assert_eq!(recs[0].alignment_start().unwrap().unwrap().get(), 51);
                assert_eq!(recs[0].reference_sequence_id(&header).unwrap().unwrap(), 0);
            }
        });

        let t2 = std::thread::spawn(move || {
            let header = aligner2.index.0.header.clone();
            for _ in 0..100000 {
                let recs = aligner2.align_read(NAME, ERCC_READ_2, ERCC_QUAL_2);
                assert_eq!(recs.len(), 1);
                assert_eq!(recs[0].alignment_start().unwrap().unwrap().get(), 501);
                assert_eq!(recs[0].reference_sequence_id(&header).unwrap().unwrap(), 0);
            }
        });

        assert!(t1.join().is_ok());
        assert!(t2.join().is_ok());
    }

    #[test]
    fn test_header() {
        let opts = StarOpts::new(ERCC_REF);
        let index = StarIndex::load(opts).unwrap();

        let mut writer = sam::io::Writer::new(Vec::new());
        writer.write_header(&index.header()).unwrap();
        let header_string = String::from_utf8(writer.get_ref().clone()).unwrap();
        assert!(header_string.starts_with("@SQ\tSN:ERCC-00002\tLN:1061\n"));
        assert!(header_string.ends_with("@SQ\tSN:ERCC-00171\tLN:505\n"));
    }
}