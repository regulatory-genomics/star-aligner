use noodles::{
    fastq,
    sam::{
        self,
        header::record::value::{map::Program, Map},
    },
};
use rayon::prelude::*;
use star_aligner::{StarAligner, StarOpts};
use std::{
    fs::File,
    io::{self, BufReader},
};

const ERCC_PARA300_PATH: &str = "/data/wenjie/test.2w.fastq";
const ERCC_REF: &str = "/data/wenjie/star-aligner/star-aligner/test/ercc92-1.2.0/star/";
fn main() {
    // get threads num from parameter
    let input = std::env::args().nth(1).unwrap();
    let threads: usize = input.parse().unwrap();
    let threads = threads as u16;

    let mut fq_reader = File::open(ERCC_PARA300_PATH)
        .map(BufReader::new)
        .map(fastq::io::Reader::new)
        .unwrap();
    let records: Vec<_> = fq_reader.records().map(Result::unwrap).collect();

    // Single end
    let opts = StarOpts::new(ERCC_REF);
    let aligner: StarAligner = StarAligner::new(opts).unwrap().with_num_threads(threads);
    let res = aligner.align_reads(&records).collect::<Vec<_>>();

    let stdout = io::stdout().lock();
    let mut writer = sam::io::Writer::new(stdout);

    let header = sam::Header::builder()
        .set_header(Default::default())
        .add_program("noodles-sam", Map::<Program>::default())
        .add_comment("an example SAM written by noodles-sam")
        .build();

    writer.write_header(&header).unwrap();
    for res_vec in res {
        for record in res_vec {
            writer.write_record(&header, &record).unwrap();
        }
    }
}
