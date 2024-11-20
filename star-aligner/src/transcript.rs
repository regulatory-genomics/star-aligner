use anyhow::{anyhow, ensure, Context, Result};
use indexmap::IndexMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::ops::Deref;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Transcriptome(Vec<Transcript>);

impl Deref for Transcriptome {
    type Target = [Transcript];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl Transcriptome {
    pub fn from_path<P: AsRef<Path>>(dir: P) -> Result<Self> {
        let genes = read_star_tab_file(dir.as_ref().join("geneInfo.tab"))?;
        let transcripts = read_star_tab_file(dir.as_ref().join("transcriptInfo.tab"))?;
        let exons = read_star_tab_file(dir.as_ref().join("exonInfo.tab"))?;
        let refseq = ReferenceSequences::new(
            dir.as_ref().join("chrNameLength.txt"),
            dir.as_ref().join("chrStart.txt"),
        )?;
        Ok(Self::new(genes, transcripts, exons, &refseq))
    }

    fn new(
        genes: Vec<StarGene>,
        transcripts: Vec<StarTranscript>,
        exons: Vec<StarExon>,
        refseq: &ReferenceSequences,
    ) -> Self {
        let new_transcripts: Vec<_> = transcripts
            .iter()
            .map(|transcript| transcript.to_transcript(&genes, &exons, refseq).unwrap())
            .collect();

        Self(new_transcripts)
    }
}

#[derive(Debug, Clone)]
pub struct Transcript {
    pub id: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64, // exclusive
    pub strand: Strand,
    pub exons: Vec<Exon>,
    pub gene_id: String,
    pub gene_name: String,
}

#[derive(Debug, Eq, PartialEq, PartialOrd, Ord, Clone)]
pub struct Exon {
    pub start: u64,
    pub end: u64, // exclusive
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

impl TryFrom<u8> for Strand {
    type Error = anyhow::Error;

    fn try_from(n: u8) -> Result<Self> {
        match n {
            0 => Ok(Self::Unknown),
            1 => Ok(Self::Forward),
            2 => Ok(Self::Reverse),
            _ => Err(anyhow!("invalid strand: {}", n)),
        }
    }
}

struct ReferenceSequences {
    sequences: IndexMap<String, u64>,
    starts: Vec<u64>,
}

impl ReferenceSequences {
    fn new<P: AsRef<Path>>(chr_name_length: P, chr_starts: P) -> Result<Self> {
        let mut buf = String::new();

        File::open(chr_name_length)?.read_to_string(&mut buf)?;
        let sequences = buf
            .lines()
            .map(|line| {
                let mut fields = line.split_whitespace();
                let name = fields
                    .next()
                    .ok_or_else(|| anyhow!("missing name"))?
                    .to_string();
                let length: u64 = fields
                    .next()
                    .ok_or_else(|| anyhow!("missing length"))?
                    .parse()?;
                Ok((name, length))
            })
            .collect::<Result<IndexMap<_, _>>>()?;

        buf.clear();
        File::open(chr_starts)?.read_to_string(&mut buf)?;
        let starts = buf
            .lines()
            .map(|line| Ok(line.parse()?))
            .collect::<Result<Vec<_>>>()?;

        ensure!(
            sequences.len() + 1 == starts.len(),
            "mismatched number of sequences and starts"
        );
        ensure!(starts.is_sorted(), "starts are not sorted");
        Ok(Self { sequences, starts })
    }
}

/// Genomic coordinates in a linear concatenated reference genome.
struct Coordinate(u64);

impl From<u64> for Coordinate {
    fn from(pos: u64) -> Self {
        Self(pos)
    }
}

impl Coordinate {
    /// Convert a coordinate to a genomic locus.
    fn to_locus<'a>(&self, refseqs: &'a ReferenceSequences) -> Locus<'a> {
        match refseqs.starts.binary_search(&self.0) {
            Ok(i) => {
                let chr = refseqs.sequences.get_index(i).unwrap().0;
                Locus {
                    chrom: chr,
                    pos: 0,
                }
            }
            Err(j) => {
                let i = j - 1;
                let chr = refseqs.sequences.get_index(i).unwrap().0;
                Locus {
                    chrom: chr,
                    pos: self.0 - refseqs.starts[i],
                }
            }
        }
    }
}

struct Locus<'a> {
    chrom: &'a str,
    pos: u64,
}

/// A STAR gene.
/// The order of these fields must match the order of the columns in the file
/// "geneInfo.tab" produced by STAR, which has no header.
struct StarGene {
    id: String,
    name: String,
    _bio_type: Option<String>,
}

impl TryFrom<&str> for StarGene {
    type Error = anyhow::Error;

    fn try_from(line: &str) -> Result<Self> {
        let mut fields = line.split('\t');
        let id = fields
            .next()
            .ok_or_else(|| anyhow!("missing id"))?
            .to_string();
        let name = fields.next().unwrap_or(&id).to_string();
        let _bio_type = fields.next().map(|s| s.to_string());
        Ok(Self { id, name, _bio_type })
    }
}

/// A STAR transcript.
/// The order of these fields must match the order of the columns in the file
/// "transcriptInfo.tab" produced by STAR, which has no header.
struct StarTranscript {
    id: String,       // transcript ID
    start: u64,       // start
    end: u64,         // end (inclusive)
    _max_end: u64,     // end max (inclusive)
    strand: u8,       // strand (0: undefined, 1: +, 2: -)
    num_exons: usize, // number of exons
    exon_idx: usize,  // index of the first exon for each transcript in exSE
    gene_idx: usize,  // transcript to gene correspondence
}

impl StarTranscript {
    fn to_transcript(&self, genes: &[StarGene], exons: &[StarExon], refseq: &ReferenceSequences) -> Result<Transcript> {
        let locus_start = Coordinate(self.start).to_locus(refseq);
        let locus_end = Coordinate(self.end).to_locus(refseq);
        let strand = Strand::try_from(self.strand)?;
        let gene_id = genes[self.gene_idx].id.clone();
        let gene_name = genes[self.gene_idx].name.clone();
        let exons = exons[self.exon_idx..self.exon_idx + self.num_exons]
            .iter()
            .map(|exon| Exon {
                start: locus_start.pos + exon.start,
                end: locus_start.pos + exon.end + 1,
            })
            .collect::<Vec<_>>();
        let expected_end = locus_end.pos + 1;
        let actual_end = exons.iter().last().unwrap().end;
        ensure!(expected_end == actual_end, "mismatched end: {} != {}", expected_end, actual_end);
        Ok(Transcript {
            id: self.id.clone(),
            chrom: locus_start.chrom.to_string(),
            start: locus_start.pos,
            end: locus_end.pos + 1,
            strand,
            exons,
            gene_id,
            gene_name,
        })
    }
}

impl TryFrom<&str> for StarTranscript {
    type Error = anyhow::Error;

    fn try_from(line: &str) -> Result<Self> {
        let mut fields = line.split('\t');
        let id = fields
            .next()
            .ok_or_else(|| anyhow!("missing id"))?
            .to_string();
        let start: u64 = fields
            .next()
            .ok_or_else(|| anyhow!("missing start"))?
            .parse()?;
        let end: u64 = fields
            .next()
            .ok_or_else(|| anyhow!("missing end"))?
            .parse()?;
        let _max_end: u64 = fields
            .next()
            .ok_or_else(|| anyhow!("missing max_end"))?
            .parse()?;
        let strand: u8 = fields
            .next()
            .ok_or_else(|| anyhow!("missing strand"))?
            .parse()?;
        let num_exons: usize = fields
            .next()
            .ok_or_else(|| anyhow!("missing num_exons"))?
            .parse()?;
        let exon_idx: usize = fields
            .next()
            .ok_or_else(|| anyhow!("missing exon_idx"))?
            .parse()?;
        let gene_idx: usize = fields
            .next()
            .ok_or_else(|| anyhow!("missing gene"))?
            .parse()?;
        Ok(Self {
            id,
            start,
            end,
            _max_end,
            strand,
            num_exons,
            exon_idx,
            gene_idx,
        })
    }
}

/// A STAR exon. The positions are 0-based and inclusive, relative to the transcript start.
#[derive(Debug)]
struct StarExon {
    start: u64,
    end: u64, // inclusive
    _cum_len: u64,
}

impl TryFrom<&str> for StarExon {
    type Error = anyhow::Error;

    fn try_from(line: &str) -> Result<Self> {
        let mut fields = line.split('\t');
        let start: u64 = fields
            .next()
            .ok_or_else(|| anyhow!("missing start"))?
            .parse()?;
        let end: u64 = fields
            .next()
            .ok_or_else(|| anyhow!("missing end"))?
            .parse()?;
        let _cum_len: u64 = fields
            .next()
            .ok_or_else(|| anyhow!("missing cum_len"))?
            .parse()?;
        Ok(Self { start, end, _cum_len })
    }
}

#[inline]
fn read_star_tab_file<P, T>(path: P) -> Result<Vec<T>>
where
    P: AsRef<Path>,
    T: for<'a> TryFrom<&'a str, Error = anyhow::Error>,
{
    let file = File::open(path.as_ref())
        .with_context(|| format!("failed to open file {:?}", path.as_ref()))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    if let Some(n) = lines.next() {
        let n: usize = n?.parse()?;
        let result = lines
            .map(|line| T::try_from(&line?))
            .collect::<Result<Vec<T>>>()?;
        ensure!(
            result.len() == n,
            "expected {} lines, got {}",
            n,
            result.len()
        );
        Ok(result)
    } else {
        Err(anyhow!("empty file"))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    const ERCC_REF: &str = "/data/Public/STAR_reference/refdata-gex-GRCm39-2024-A/star_2.7.11b";

    #[test]
    fn test_star() {
        let transcriptome = Transcriptome::from_path(ERCC_REF).unwrap();
        println!("{:?}", transcriptome[999]);
    }
}