use rayon::prelude::*;

use crate::args::ScoreArgs;
use crate::common::{encode_dna, open_file_or_stdout};
use crate::constants::*;
use crate::errors::*;
use crate::index::KMerIndex;
use crate::progress;
use crate::score::calculate_score;
use crate::table;

fn build_row(index: &KMerIndex, idx: usize, row: &mut Vec<String>) {
    let pam = &index.enzyme().pam;
    let value = row
        .first()
        .expect("unexpected empty table row")
        .to_ascii_uppercase();

    if value.len() >= KMER_LEN + pam.len() && pam.matches(value.as_bytes()) {
        let start = value.len() - pam.len() - KMER_LEN;
        let kmer = &value.as_bytes()[start..start + KMER_LEN];
        assert!(kmer.len() == KMER_LEN);

        if let Some(kmer) = encode_dna(kmer) {
            row.push(calculate_score(&index, kmer).to_string());
            return;
        }
    }

    // Not a valid gRNA sequence; either a header or (presumably) DNA containing Ns
    if idx == 0 {
        row.push("Score".into());
    } else {
        row.push("NA".into());
    }
}

pub fn main(args: &ScoreArgs) -> Result<()> {
    ::rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .chain_err(|| "failed to build thread pool")?;

    eprintln!("\nReading K-mers from {:?}", args.index);
    let index = KMerIndex::read(&args.index)
        .chain_err(|| format!("failed to read K-mer index {:?}", args.index))?;
    eprintln!("  {}", index.summarize());

    eprintln!("Reading target sites from {:?}", args.table);
    let mut table =
        table::read(&args.table).chain_err(|| "failed to read table of target sites")?;
    eprintln!("  read {} target sites from table.", table.len());

    let progress = progress::default(table.len());
    table.par_iter_mut().enumerate().for_each(|(idx, row)| {
        build_row(&index, idx, row);

        progress.inc(1);
    });

    progress.finish();

    let mut out = open_file_or_stdout(&args.output)?;
    for row in table {
        writeln!(out, "{}", row.join("\t")).chain_err(|| "failed to write output row")?;
    }

    Ok(())
}
