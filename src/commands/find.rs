use bio::alphabets::dna;
use bio::io::fasta::Reader;
use bio_types::strand::Strand;
use rayon::prelude::*;

use crate::args::FindArgs;
use crate::common::open_file_or_stdout;
use crate::errors::*;
use crate::index::KMerIndex;
use crate::pam::Position;
use crate::progress;
use crate::progress::ProgressBar;
use crate::score::calculate_score;

#[derive(Debug)]
struct TargetSite {
    start: isize,
    end: isize,
    cutsite: Option<isize>,
    strand: Strand,
    sequence: Vec<u8>,
    score: u64,
    depth: usize,
}

type TargetSites = Vec<TargetSite>;

fn collect_forward_targets(sequence: &[u8], index: &KMerIndex, pg: &ProgressBar) -> TargetSites {
    let enzyme = index.enzyme();
    let pam = &enzyme.pam;
    let cutsite = enzyme.cutsite;

    sequence
        .par_windows(enzyme.grna_len)
        .enumerate()
        .filter_map(|(idx, window)| {
            let result = if let Some((pam_pos, kmer)) = pam.kmer(window) {
                Some(TargetSite {
                    start: idx as isize,
                    end: idx as isize + window.len() as isize,
                    cutsite: cutsite.map(|v| (idx + pam_pos) as isize + v as isize),
                    strand: Strand::Forward,
                    sequence: window.to_owned(),
                    score: calculate_score(index, kmer),
                    depth: 0,
                })
            } else {
                None
            };

            pg.inc(1);
            result
        })
        .collect()
}

fn collect_reverse_targets(sequence: &[u8], index: &KMerIndex, pg: &ProgressBar) -> TargetSites {
    let sequence = dna::revcomp(sequence);
    let mut sites = collect_forward_targets(&sequence, index, pg);

    for site in &mut sites {
        let start = sequence.len() as isize - site.end;
        let end = sequence.len() as isize - site.start;

        site.start = start;
        site.end = end;
        site.cutsite = site.cutsite.map(|v| sequence.len() as isize - v);
        site.strand = Strand::Reverse;
    }

    sites
}

fn collect_targets(index: &KMerIndex, sequence: &[u8], pg: &ProgressBar) -> TargetSites {
    let mut targets = Vec::new();
    targets.append(&mut collect_forward_targets(&sequence, index, pg));
    targets.append(&mut collect_reverse_targets(&sequence, index, pg));
    targets.sort_unstable_by_key(|v| v.cutsite);

    targets
}

pub fn main(args: &FindArgs) -> Result<()> {
    ::rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .chain_err(|| "failed to build thread pool")?;

    eprintln!("\nReading K-mers from {:?}", args.index);
    let index = KMerIndex::read(&args.index)
        .chain_err(|| format!("failed to read K-mer index {:?}", &args.index))?;
    eprintln!("  {}", index.summarize());

    eprintln!("Finding target sites in {:?}", &args.targets);
    let reader = Reader::from_file(&args.targets)
        .chain_err(|| format!("failed to open FASTA file {:?}", args.targets))?;

    let mut out = open_file_or_stdout(&args.output)?;
    writeln!(out, "Sequence\tName\tStart\tEnd\tCutsite\tStrand\tScore")
        .chain_err(|| "failed to write output header")?;

    let enzyme = index.enzyme();
    let grna_len = enzyme.grna_len;
    let pam_len = enzyme.pam.len();
    let pam_offset = match enzyme.pam.position() {
        Position::Head => 0,
        Position::Tail => enzyme.grna_len - pam_len,
    };

    for (idx, target) in reader.records().enumerate() {
        let target = target.chain_err(|| "failed to read sequence")?;
        target.check().map_err(|v| ErrorKind::Msg(v.into()))?;

        let sequence = target.seq().to_ascii_uppercase();
        let num_sites = usize::max(sequence.len(), grna_len - 1) - grna_len - 1;
        let prefix = format!("{}. {}: ", idx + 1, target.id());
        let progress = progress::with_prefix(num_sites * 2, &prefix);

        for site in collect_targets(&index, &sequence, &progress) {
            let mut sequence = site.sequence.clone();
            for nuc in &mut sequence[pam_offset..pam_offset + pam_len] {
                *nuc = nuc.to_ascii_lowercase();
            }

            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                String::from_utf8_lossy(&sequence),
                target.id(),
                site.start + 1,
                site.end,
                site.cutsite
                    .map(|v| (v + 1).to_string())
                    .unwrap_or_else(|| "NA".to_owned()),
                site.strand.strand_symbol(),
                site.score,
            )
            .chain_err(|| "failed to write output row")?;
        }

        progress.finish();
    }

    Ok(())
}
