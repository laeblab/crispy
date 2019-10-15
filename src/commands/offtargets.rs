use std::collections::HashSet;

use crate::args::OffTargetsArgs;
use crate::common::encode_dna;
use crate::constants::*;
use crate::errors::*;
use crate::index::KMerIndex;
use crate::pam::Position;
use crate::score::find_offtargets;
use crate::table;

fn print_off_targets(index: &KMerIndex, query: &str, value: &str) {
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

            for position in find_offtargets(index, kmer) {
                let (offset_start, offset_end) = if position.strand() == '+' {
                    (offset_start, offset_end)
                } else {
                    (1 - offset_end, 1 - offset_start)
                };

                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}",
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
                );
            }

            return;
        }
    }

    eprint!("WARNING: Could not look up off-targets for '{:?}'", value);
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

    println!("Query\tName\tStart\tEnd\tCutsite\tStrand");

    let mut values = HashSet::new();
    for input_row in &table {
        let value = input_row.first().expect("unexpected empty table row");
        let query = value.to_ascii_uppercase();

        if values.insert(value.clone()) {
            print_off_targets(&index, &query, &value);
        }
    }

    Ok(())
}
