use std::fmt::Debug;
use std::io::prelude::*;
use std::path::Path;

use bio::alphabets::dna;
use bio::io::fasta::Reader;

use crate::args::IndexArgs;
use crate::constants::*;
use crate::enzyme::Enzyme;
use crate::errors::*;
use crate::index::{KMerIndex, KMerMap, Position};

fn collect_hashes<P: AsRef<Path> + Debug>(
    filename: &P,
    enzyme: &Enzyme,
    positions: bool,
) -> Result<(Vec<String>, KMerMap)> {
    eprintln!("Finding target sites in {:?}", &filename);
    let file = Reader::from_file(&filename)
        .chain_err(|| format!("failed to open FASTA file {:?}", filename))?;

    let mut refseqs = Vec::new();
    let mut hashes = if positions {
        KMerMap::Positions(vec![Vec::new(); KMER_COUNT])
    } else {
        KMerMap::Counts(vec![0; KMER_COUNT])
    };

    let mut running_size = 0;
    let timer = ::std::time::Instant::now();
    for (refseq, record) in file.records().enumerate() {
        let record = record.chain_err(|| "failed to read FASTA sequence")?;
        record.check().map_err(|v| ErrorKind::Msg(v.into()))?;

        refseqs.push(record.id().to_owned());

        let sequence = record.seq().to_ascii_uppercase();
        for (idx, window) in sequence.windows(enzyme.pam.len() + KMER_LEN).enumerate() {
            if let Some((pam_pos, kmer)) = enzyme.pam.kmer(window) {
                // If the cut-site is unknown, then save the starting position of the target seq
                let pos = ((idx + pam_pos) as isize + enzyme.cutsite.unwrap_or(0)) as i32;

                hashes.add(Position::forward(refseq as u32, pos), kmer);
            }
        }

        for (idx, window) in dna::revcomp(&sequence)
            .windows(enzyme.pam.len() + KMER_LEN)
            .enumerate()
        {
            if let Some((pam_pos, kmer)) = enzyme.pam.kmer(window) {
                // If the cut-site is unknown, then save the starting position of the target seq
                let pos = sequence.len() as i32
                    - ((idx + pam_pos) as isize + enzyme.cutsite.unwrap_or(0)) as i32;

                hashes.add(Position::reverse(refseq as u32, pos), kmer);
            }
        }

        let seconds = timer.elapsed().as_secs() as usize;
        running_size += sequence.len();

        print!(
            "\r  Processed {} Mbp in {} seconds ({:.1} Mbp/s)",
            running_size / 1_000_000,
            seconds,
            (running_size / 1_000_000) as f64 / ::std::cmp::max(1, seconds) as f64
        );
        ::std::io::stdout().flush().expect("unable to flush STDOUT");
    }

    println!();

    Ok((refseqs, hashes))
}

pub fn main(args: &IndexArgs) -> Result<()> {
    println!("Finding targets in {:?}", &args.fasta);
    let (refseqs, hashes) = collect_hashes(&args.fasta, &args.enzyme, args.positions)
        .chain_err(|| "failed to collect target sequence frequencies")?;

    let filename = match &args.output {
        Some(filename) => filename.clone(),
        None => args.fasta.clone() + args.enzyme.extension,
    };

    let index = KMerIndex::new(&args.enzyme, refseqs, hashes);
    println!("  Unique targets found: {}", index.kmer_count());
    println!("  Writing list of target frequencies to {:?}", filename);
    index
        .write(&filename)
        .chain_err(|| "failed to write target site frequencies")?;

    Ok(())
}
