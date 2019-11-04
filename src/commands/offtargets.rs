use std::collections::HashSet;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use bio::alphabets::dna;
use bio::io::fasta::IndexedReader;

use crate::args::OffTargetsArgs;
use crate::common::{encode_dna, open_file_or_stdout};
use crate::constants::*;
use crate::enzyme::Enzyme;
use crate::errors::*;
use crate::index::KMerIndex;
use crate::pam::Position;
use crate::score::find_offtargets;
use crate::table;

struct OfftargetReader {
    reader: Option<IndexedReader<File>>,
}

impl OfftargetReader {
    fn new(index: &str, fasta: &Option<String>) -> Result<Self> {
        let fasta_path = match fasta {
            Some(path) => path,
            None => &index[..index.rfind('.').unwrap_or(index.len())],
        };

        let reader = if Path::new(fasta_path).exists() {
            let mut fai_path = fasta_path.to_owned();
            fai_path.push_str(".fai");

            if Path::new(&fai_path).exists() {
                IndexedReader::from_file(&fasta_path)
                    .map(|v| Some(v))
                    .chain_err(|| format!("failed to read FASTA file {:?}", fasta_path))?
            } else {
                eprintln!("{:?} not indexed; cannot fetch off-target seqs", fasta_path);
                eprintln!("Please run: samtools faidx '{}'", fasta_path);

                None
            }
        } else {
            eprintln!("FASTA file not found at {:?}", fasta_path);

            None
        };

        Ok(OfftargetReader { reader })
    }

    fn fetch(
        &mut self,
        refseq: &str,
        mut start: isize,
        end: isize,
        strand: char,
    ) -> Result<Option<Vec<u8>>> {
        assert!(start < end);

        if let Some(reader) = &mut self.reader {
            let len = (end - start) as usize;

            if end <= 0 {
                Ok(Some(vec![b'N'; len]))
            } else {
                let mut seq = Vec::with_capacity(len);

                // 3' PAM gRNAs may extend past the beginning of the refseq
                if start < 0 {
                    seq.resize(start.abs() as usize, b'N');
                    start = 0;
                }

                reader
                    .fetch(refseq, start as u64, end as u64)
                    .chain_err(|| "failed to fetch refseq")?;
                reader
                    .read(&mut seq)
                    .chain_err(|| "failed to read refseq")?;

                // 5' PAM gRNAs may extend past the end of the refseq
                seq.resize(len, b'N');

                if strand == '-' {
                    seq = dna::revcomp(seq);
                }

                Ok(Some(seq))
            }
        } else {
            Ok(None)
        }
    }
}

/// Indicates the KMer in uppercase, rendering everything else in lowercase
fn format_guide_rna(enzyme: &Enzyme, seq: &[u8]) -> String {
    let mut seq = seq.to_ascii_lowercase();
    for c in enzyme.pam.kmer_slice_mut(&mut seq) {
        *c = c.to_ascii_uppercase();
    }

    String::from_utf8_lossy(&seq).to_string()
}

fn write_off_targets(
    out: &mut dyn Write,
    fasta: &mut OfftargetReader,
    index: &KMerIndex,
    query: &str,
    value: &str,
    min_score: u64,
) -> Result<()> {
    let enzyme = index.enzyme();
    let refseqs = index.refseqs();
    let pam = &enzyme.pam;

    if query.len() >= KMER_LEN + pam.len() && pam.matches(query.as_bytes()) {
        let query = query.as_bytes();
        let kmer = enzyme.pam.kmer_slice(query).to_ascii_uppercase();
        assert!(kmer.len() == KMER_LEN);

        if let Some(kmer) = encode_dna(&kmer) {
            let grna_len = enzyme.grna_len as isize;
            let pam_len = pam.len() as isize;
            let cutsite = enzyme.cutsite.unwrap_or(0);

            let (offset_start, offset_end) = match pam.position() {
                Position::Head => (1 - cutsite, grna_len - cutsite),
                Position::Tail => (-cutsite - grna_len + pam_len + 1, pam_len - cutsite),
            };

            for (score, position) in find_offtargets(index, kmer, min_score) {
                let (offset_start, offset_end) = if position.strand() == '+' {
                    (offset_start, offset_end)
                } else {
                    (1 - offset_end, 1 - offset_start)
                };

                let refseq = &refseqs[position.refseq() as usize];
                let start = position.pos() as isize + offset_start;
                let end = position.pos() as isize + offset_end;
                // fetch uses 0-based start, 1-based end
                let offtarget = fasta
                    .fetch(refseq, start - 1, end, position.strand())
                    .chain_err(|| "failed to fetch offtarget sequence")?;

                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    format_guide_rna(index.enzyme(), query),
                    if let Some(seq) = offtarget {
                        format_guide_rna(index.enzyme(), &seq)
                    } else {
                        "NA".to_owned()
                    },
                    refseq,
                    start,
                    end,
                    if enzyme.cutsite.is_none() {
                        "NA".to_owned()
                    } else {
                        (position.pos() + 1).to_string()
                    },
                    position.strand(),
                    score,
                )
                .chain_err(|| "failed to write output row")?;
            }

            return Ok(());
        }
    }

    eprintln!("WARNING: Could not look up off-targets for '{:?}'", value);

    Ok(())
}

pub fn main(args: &OffTargetsArgs) -> Result<()> {
    eprintln!("\nReading K-mers from {:?}", args.index);
    let index = KMerIndex::read(&args.index)
        .chain_err(|| format!("failed to read K-mer index {:?}", args.index))?;
    eprintln!("  {}", index.summarize());

    if !index.has_positions() {
        return Err("fasta not indexed with --positions; cannot find off-targets ".into());
    }

    eprintln!("Reading target sites from {:?}", args.table);
    let table = table::read(&args.table).chain_err(|| "failed to read table of target sites")?;
    eprintln!("  read {} target sites from table.", table.len());

    let mut reader =
        OfftargetReader::new(&args.index, &args.fasta).chain_err(|| "failed to open FASTA file")?;

    let mut out = open_file_or_stdout(&args.output)?;
    writeln!(
        out,
        "Query\tOfftarget\tName\tStart\tEnd\tCutsite\tStrand\tScore"
    )
    .chain_err(|| "failed to write output header")?;

    let mut values = HashSet::new();
    for input_row in &table {
        let value = input_row.first().expect("unexpected empty table row");
        let query = value.to_ascii_uppercase();

        if values.insert(value.clone()) {
            write_off_targets(
                &mut out,
                &mut reader,
                &index,
                &query,
                &value,
                args.min_score,
            )?;
        }
    }

    Ok(())
}
