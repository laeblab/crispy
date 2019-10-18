use std::collections::HashSet;
use std::io::Write;

use crate::args::OffTargetsArgs;
use crate::common::{encode_dna, open_file_or_stdout};
use crate::constants::*;
use crate::errors::*;
use crate::index::KMerIndex;
use crate::pam::Position;
use crate::score::find_offtargets;
use crate::table;

fn write_off_targets(
    out: &mut dyn Write,
    index: &KMerIndex,
    query: &str,
    value: &str,
    min_score: u64,
) -> Result<()> {
    let enzyme = index.enzyme();
    let refseqs = index.refseqs();
    let pam = &enzyme.pam;

    if query.len() >= KMER_LEN + pam.len() && pam.matches(query.as_bytes()) {
        let start = match pam.position() {
            Position::Head => pam.len(),
            Position::Tail => query.len() - pam.len() - KMER_LEN,
        };

        let kmer = query.as_bytes()[start..start + KMER_LEN].to_ascii_uppercase();
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

                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    query,
                    refseqs[position.refseq() as usize],
                    position.pos() as isize + offset_start,
                    position.pos() as isize + offset_end,
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

    eprint!("WARNING: Could not look up off-targets for '{:?}'", value);

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

    let mut out = open_file_or_stdout(&args.output)?;
    writeln!(out, "Query\tName\tStart\tEnd\tCutsite\tStrand\tScore")
        .chain_err(|| "failed to write output header")?;

    let mut values = HashSet::new();
    for input_row in &table {
        let value = input_row.first().expect("unexpected empty table row");
        let query = value.to_ascii_uppercase();

        if values.insert(value.clone()) {
            write_off_targets(&mut out, &index, &query, &value, args.min_score)?;
        }
    }

    Ok(())
}
