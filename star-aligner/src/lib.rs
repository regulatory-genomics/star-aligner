use anyhow::{bail, Result};
use noodles::fastq::Record as FastqRecord;
use noodles::sam::{
    header::record::value::{map::ReferenceSequence, Map},
    Header, Record as SamRecord,
};
use noodles::{fastq, sam};
use rayon::prelude::*;
use std::sync::Mutex;
use std::{
    // collections::HashMap,
    ffi::{c_char, c_int, CStr, CString},
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    sync::Arc,
    vec,
};

pub trait TranscriptomeAligner {
    fn chunk_size(&self) -> usize;

    fn header(&self) -> sam::Header;

    // diff with genome aligner, return a vector of sam records
    fn align_reads(
        &mut self,
        records: &mut [fastq::Record],
    ) -> impl ExactSizeIterator<Item = Vec<sam::Record>>;

    fn align_read_pairs(
        &mut self,
        records: &mut [(fastq::Record, fastq::Record)],
    ) -> impl ExactSizeIterator<Item = Vec<(sam::Record, sam::Record)>>;
}

struct InnerStarIndex {
    index: *const star_sys::StarRef,
    opts: StarOpts,
    header: Header,
}

// Ensure InnerStarIndex is Send and Sync
unsafe impl Send for InnerStarIndex {}
unsafe impl Sync for InnerStarIndex {}

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

    pub fn get_aligner(&self) -> StarAligner {
        StarAligner::new(self.clone())
    }
}

/// StarOpts contains the parameters which will be used for the STAR aligner.
/// Currently the array of argument strings is passed directly
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

/// StarAligner aligns single reads or read-pairs to the reference it is initialized with, and returns
/// rust_htslib Record objects
pub struct StarAligner {
    aligner: *mut star_sys::Aligner,
    index: StarIndex,
    // sam_buf: SamRecord,
    // fastq1: FastqRecord,
    // fastq2: FastqRecord,
}

impl Drop for StarAligner {
    fn drop(&mut self) {
        unsafe { star_sys::destroy_aligner(self.aligner) };
    }
}

unsafe impl Send for StarAligner {}
unsafe impl Sync for StarAligner {}

impl Clone for StarAligner {
    fn clone(&self) -> StarAligner {
        StarAligner::new(self.index.clone())
    }
}

impl StarAligner {
    pub fn new(index: StarIndex) -> Self {
        let aligner = unsafe { star_sys::init_aligner_from_ref(index.0.as_ref().index) };
        Self {
            aligner,
            index,
            // sam_buf: SamRecord::default(),
            // fastq1: FastqRecord::default(),
            // fastq2: FastqRecord::default(),
        }
    }

    pub fn get_header(&self) -> &Header {
        &self.index.0.header
    }

    // align single read in parallel
    pub fn align_reads(
        &mut self,
        fqrecs: &mut [fastq::Record],
        n_threads: usize,
    ) -> Result<Vec<SamRecord>> {
        // set thread pool
        rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build_global()?;
        // split the reads into chunks by the number of threads
        let chunk_size = (fqrecs.len() + n_threads - 1) / n_threads;
        let fqrecs_chunks: Vec<_> = fqrecs.chunks(chunk_size).collect();

        // construct an Arc mut clone to share aligner among threads
        let aligner = Arc::new(Mutex::new(self.clone()));

        // align reads in parallel
        let records: Vec<SamRecord> = fqrecs_chunks
            .par_iter()
            .flat_map(|chunk| {
                chunk
                    .iter()
                    .flat_map(|fqrec| {
                        // lock the aligner
                        let mut aligner = aligner.lock().unwrap();
                        aligner.align_read(fqrec).unwrap()
                    })
                    .collect::<Vec<SamRecord>>()
            })
            .collect();

        Ok(records)
    }

    // align read pairs in parallel
    pub fn align_pair_reads(
        &mut self,
        fqrecs: &mut [(fastq::Record, fastq::Record)],
        n_threads: usize,
    ) -> Result<Vec<(SamRecord, SamRecord)>> {
        // set thread pool
        rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build_global()?;
        // split the reads into chunks by the number of threads
        let chunk_size = (fqrecs.len() + n_threads - 1) / n_threads;
        let fqrecs_chunks: Vec<_> = fqrecs.chunks(chunk_size).collect();

        // construct an Arc mut clone to share aligner among threads
        let aligner = Arc::new(Mutex::new(self.clone()));

        // align reads in parallel
        let records: Vec<(SamRecord, SamRecord)> = fqrecs_chunks
            .par_iter()
            .flat_map(|chunk| {
                chunk
                    .iter()
                    .map(|pairfqrec| {
                        // lock the aligner
                        let mut aligner = aligner.lock().unwrap();
                        aligner
                            .align_read_pair((&pairfqrec.0, &pairfqrec.1))
                            .unwrap()
                    })
                    .flat_map(|(first, second)| first.into_iter().zip(second))
                    .collect::<Vec<(SamRecord, SamRecord)>>()
            })
            .collect();

        Ok(records)
    }

    pub fn align_read(&mut self, fqrec: &FastqRecord) -> Result<Vec<SamRecord>> {
        // init a aligner
        self.aligner = unsafe { star_sys::init_aligner_from_ref(self.index.0.as_ref().index) };

        let mut fq_buf = Vec::new();
        fqrec_to_buf(&mut fq_buf, fqrec)?;
        let fq_cstr = CStr::from_bytes_with_nul(&fq_buf)?;
        let res: *const c_char = unsafe { star_sys::align_read(self.aligner, fq_cstr.as_ptr()) };
        if res.is_null() {
            bail!("STAR returned null alignment");
        }

        let fqname = fqrec.name();

        let cstr = unsafe { CStr::from_ptr(res) };
        // aln_buf.clear();
        // aln_buf.extend_from_slice(cstr.to_bytes());
        let aln_buf = cstr.to_bytes();

        let mut sam_buf = Vec::new();

        let mut records = Vec::new();
        for slc in aln_buf.split(|c| *c == b'\n') {
            if !slc.is_empty() {
                sam_buf.clear();
                sam_buf.extend_from_slice(fqname);
                sam_buf.extend_from_slice(slc);
                let record = sam_buf.as_slice().try_into().unwrap();
                records.push(record);
            }
        }

        unsafe {
            libc::free(res as *mut libc::c_void);
        }
        Ok(records)
    }

    pub fn align_read_pair(
        &mut self,
        fqrec_pair: (&FastqRecord, &FastqRecord),
    ) -> Result<(Vec<SamRecord>, Vec<SamRecord>)> {
        let (fqrec1, fqrec2) = fqrec_pair;

        // convert FastqRecord to CStr
        let mut fq1_buf = Vec::new();
        let mut fq2_buf = Vec::new();
        fqrec_to_buf(&mut fq1_buf, fqrec1)?;
        fqrec_to_buf(&mut fq2_buf, fqrec2)?;
        let fq1cstr = CStr::from_bytes_with_nul(&fq1_buf)?;
        let fq2cstr = CStr::from_bytes_with_nul(&fq2_buf)?;

        let res: *const c_char =
            unsafe { star_sys::align_read_pair(self.aligner, fq1cstr.as_ptr(), fq2cstr.as_ptr()) };
        if res.is_null() {
            bail!("STAR returned null alignment");
        }

        // let fqname2 = fqrec2.name();
        let fqname = fqrec1.name();

        let rescstr = unsafe { CStr::from_ptr(res) };
        let aln_buf = rescstr.to_bytes();
        let mut sam_buf = Vec::new();

        let mut full_records: Vec<SamRecord> = Vec::new();
        for slc in aln_buf.split(|c| *c == b'\n') {
            if !slc.is_empty() {
                sam_buf.clear();
                sam_buf.extend_from_slice(fqname);
                sam_buf.extend_from_slice(slc);
                let record = sam_buf.as_slice().try_into().unwrap();
                full_records.push(record);
            }
        }

        unsafe {
            libc::free(res as *mut libc::c_void);
        }
        // Partition the records into first mate and second mate
        let mut first_vec: Vec<sam::Record> = Vec::new();
        let mut second_vec: Vec<sam::Record> = Vec::new();
        for rec in full_records {
            if rec.flags().unwrap().is_first_segment() {
                first_vec.push(rec);
            } else {
                second_vec.push(rec);
            }
        }
        Ok((first_vec, second_vec))
    }
}

/// Convert FastqRecord to CStr
fn fqrec_to_buf(buf: &mut Vec<u8>, fqrec: &FastqRecord) -> Result<()> {
    buf.clear();
    buf.push(b'@');
    buf.extend_from_slice(fqrec.name());
    buf.push(b'\n');
    buf.extend_from_slice(fqrec.sequence());
    buf.push(b'\n');
    buf.push(b'+');
    buf.push(b'\n');
    buf.extend_from_slice(fqrec.quality_scores());
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

    use fastq::record::Definition;

    use super::*;

    /// References to some commonly used reference genomes for testing purposes
    // pub const DEFAULT_REF_1: &str = "/mnt/opt/refdata_cellranger/mm10-3.0.0/star";
    // pub const DEFAULT_REF_2: &str = "/mnt/opt/refdata_cellranger/GRCh38-3.0.0/star";
    // pub const DEFAULT_REF_3: &str = "/mnt/opt/refdata_cellranger/GRCh38-1.2.0/star";
    // pub const DEFAULT_REF_4: &str = "//mnt/opt/refdata_cellranger/hg19_and_mm10-3.0.0/star";

    const ERCC_REF: &str = "/data/kzhang/dev/orbit/star-aligner/test/ercc92-1.2.0/star/";

    const NAME: &[u8] = b"NAME";
    const ERCC_READ_1: &[u8] = b"GCATCCAGACCGTCGGCTGATCGTGGTTTTACTAGGCTAGACTAGCGTACGAGCACTATGGTCAGTAATTCCTGGAGGAATAGGTACCAAGAAAAAAACG";
    const ERCC_QUAL_1: &[u8] = b"????????????????????????????????????????????????????????????????????????????????????????????????????";

    const ERCC_READ_2: &[u8] = b"GGAGACGAATTGCCAGAATTATTAACTGCGCAGTTAGGGCAGCGTCTGAGGAAGTTTGCTGCGGTTTCGCCTTGACCGCGGGAAGGAGACATAACGATAG";
    const ERCC_QUAL_2: &[u8] = b"????????????????????????????????????????????????????????????????????????????????????????????????????";

    const ERCC_READ_3: &[u8] = b"AACTTAATGGACGGG";
    const ERCC_QUAL_3: &[u8] = b"???????????????";

    const ERCC_READ_4: &[u8] = b"AATCCACTCAATAAATCTAAAAAC";
    const ERCC_QUAL_4: &[u8] = b"????????????????????????";

    const ERCC_PE_1_PATH: &str = "/data/wenjie/sim.R1.fq";
    const ERCC_PE_2_PATH: &str = "/data/wenjie/sim.R2.fq";

    const ERCC_PARA300_PATH: &str = "/data/wenjie/para1k.fq";
    // fn have_refs() -> bool {
    //     Path::new("/mnt/opt/refdata_cellranger").exists()
    // }

    fn fake_fqrec(name: &[u8], seq: &[u8], qual: &[u8]) -> FastqRecord {
        FastqRecord::new(Definition::new(name, ""), seq, qual)
    }

    #[test]
    fn test_old_multithreads() {
        let opts = StarOpts::new(ERCC_REF);
        let index = StarIndex::load(opts).unwrap();
        let aligner = index.get_aligner();

        let mut fqrdr = File::open(ERCC_PARA300_PATH)
            .map(BufReader::new)
            .map(fastq::io::Reader::new)
            .unwrap();
        let mut fqrecs = Vec::new();
        for fqrec in fqrdr.records() {
            let fqrec = fqrec.unwrap();
            fqrecs.push(fqrec);
        }

        let threads = 64;
        // set number of threads
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads as usize)
            .build_global()
            .unwrap();
        // let fqrecs = vec![fake_rec1, fake_rec2, fake_rec3, fake_rec4];

        let _records: Vec<Vec<SamRecord>> = fqrecs
            .par_iter()
            .map(|fqrec| {
                let mut aligner = aligner.clone();
                aligner.align_read(fqrec).unwrap()
            })
            .collect();
    }

    #[test]
    fn test_multithreads_align() {
        let opts = StarOpts::new(ERCC_REF);
        let index = StarIndex::load(opts).unwrap();
        let mut aligner = index.get_aligner();

        let mut fqrdr = File::open(ERCC_PARA300_PATH)
            .map(BufReader::new)
            .map(fastq::io::Reader::new)
            .unwrap();
        let mut fqrecs = Vec::new();
        for fqrec in fqrdr.records() {
            let fqrec = fqrec.unwrap();
            fqrecs.push(fqrec);
        }
        let _res = aligner.align_reads(&mut fqrecs, 64);
    }

    #[test]
    fn test_singlethread_align() {
        let opts = StarOpts::new(ERCC_REF);
        let index = StarIndex::load(opts).unwrap();
        let aligner = index.get_aligner();

        let mut fqrdr = File::open(ERCC_PARA300_PATH)
            .map(BufReader::new)
            .map(fastq::io::Reader::new)
            .unwrap();
        let mut fqrecs = Vec::new();
        for fqrec in fqrdr.records() {
            let fqrec = fqrec.unwrap();
            fqrecs.push(fqrec);
        }

        let _records: Vec<Vec<SamRecord>> = fqrecs
            .iter()
            .map(|fqrec| {
                let mut aligner = aligner.clone();
                aligner.align_read(fqrec).unwrap()
            })
            .collect();
    }

    #[test]
    fn test_empty_tiny_reads() {
        let opts = StarOpts::new(ERCC_REF);
        let index = StarIndex::load(opts).unwrap();
        let mut aligner = StarAligner::new(index);
        // let mut alinger = aligner.geta

        let fake_rec1 = fake_fqrec(b"fk1", b"A", b"?");

        let recs = aligner.align_read(&fake_rec1).unwrap();
        println!("{:?}", recs);
        let fake_pe_rec1 = fake_fqrec(b"fk", b"A", b"?");
        let fake_pe_rec2 = fake_fqrec(b"fk", b"C", b"?");
        let (recs1, recs2) = aligner
            .align_read_pair((&fake_pe_rec1, &fake_pe_rec2))
            .unwrap();
        println!("{:?}, {:?}", recs1, recs2);
    }

    #[test]
    fn test_ercc_align() {
        let opts = StarOpts::new(ERCC_REF);
        let index = StarIndex::load(opts).unwrap();
        let mut aligner = StarAligner::new(index);
        let header = aligner.index.0.header.clone();

        let fake_rec1 = fake_fqrec(NAME, ERCC_READ_1, ERCC_QUAL_1);

        let recs = aligner.align_read(&fake_rec1).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].alignment_start().unwrap().unwrap().get(), 51);
        assert_eq!(recs[0].reference_sequence_id(&header).unwrap().unwrap(), 0);
        println!("{:?}", recs);

        let fake_rec2 = fake_fqrec(NAME, ERCC_READ_2, ERCC_QUAL_2);
        let recs = aligner.align_read(&fake_rec2).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].reference_sequence_id(&header).unwrap().unwrap(), 0);
        assert_eq!(recs[0].alignment_start().unwrap().unwrap().get(), 501);
        println!("{:?}", recs);

        let fake_rec3 = fake_fqrec(NAME, ERCC_READ_3, ERCC_QUAL_3);
        let recs = aligner.align_read(&fake_rec3).unwrap();
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

        let fake_rec4 = fake_fqrec(NAME, ERCC_READ_4, ERCC_QUAL_4);
        let recs = aligner.align_read(&fake_rec4).unwrap();
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
    fn test_ercc_pefile_align() {
        let opts = StarOpts::new(ERCC_REF);
        let index = StarIndex::load(opts).unwrap();
        let mut aligner = StarAligner::new(index);
        // let header = aligner.index.0.header.clone();

        let mut fq1rdr = File::open(ERCC_PE_1_PATH)
            .map(BufReader::new)
            .map(fastq::io::Reader::new)
            .unwrap();
        let mut record1 = FastqRecord::default();
        let _ = fastq::Reader::read_record(&mut fq1rdr, &mut record1);
        println!("rec1:{:?}", record1);

        let mut fq2rdr = File::open(ERCC_PE_2_PATH)
            .map(BufReader::new)
            .map(fastq::io::Reader::new)
            .unwrap();
        let mut record2 = FastqRecord::default();
        let _ = fastq::Reader::read_record(&mut fq2rdr, &mut record2);
        println!("rec2:{:?}", record2);

        let (pair_res1, pair_res2) = aligner.align_read_pair((&record1, &record2)).unwrap();
        println!("{:?}", pair_res1);
        println!("{:?}", pair_res2);
    }

    #[test]
    fn test_transcriptome_min_score() {
        let mut opts = StarOpts::new(ERCC_REF);
        opts.out_filter_score_min = 20;
        let index = StarIndex::load(opts).unwrap();
        let mut aligner = StarAligner::new(index);
        let header = aligner.index.0.header.clone();

        let fake_rec = fake_fqrec(NAME, ERCC_READ_3, ERCC_QUAL_3);

        let recs = aligner.align_read(&fake_rec).unwrap();
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
                let fake_rec1 = fake_fqrec(NAME, ERCC_READ_1, ERCC_QUAL_1);
                let recs = aligner1.align_read(&fake_rec1).unwrap();
                assert_eq!(recs.len(), 1);
                assert_eq!(recs[0].alignment_start().unwrap().unwrap().get(), 51);
                assert_eq!(recs[0].reference_sequence_id(&header).unwrap().unwrap(), 0);
            }
        });

        let t2 = std::thread::spawn(move || {
            let header = aligner2.index.0.header.clone();
            for _ in 0..100000 {
                let fake_rec2 = fake_fqrec(NAME, ERCC_READ_2, ERCC_QUAL_2);
                let recs = aligner2.align_read(&fake_rec2).unwrap();
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
        writer.write_header(index.header()).unwrap();
        let header_string = String::from_utf8(writer.get_ref().clone()).unwrap();
        assert!(header_string.starts_with("@SQ\tSN:ERCC-00002\tLN:1061\n"));
        assert!(header_string.ends_with("@SQ\tSN:ERCC-00171\tLN:505\n"));
    }
}
