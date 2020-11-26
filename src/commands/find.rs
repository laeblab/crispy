use bio::alphabets::dna;
use bio::io::bed;
use bio::io::fasta;
use bio_types::strand::Strand;
use rayon::prelude::*;
use std::collections::HashMap;

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
    cutsite: isize,
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
                    cutsite: (idx + pam_pos) as isize + cutsite,
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
        site.cutsite = sequence.len() as isize - site.cutsite;
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

fn print_target(
    record: &bed::Record,
    site: &TargetSite,
    sequence: &[u8],
    offset: isize,
    out: &mut Box<dyn std::io::Write>,
) -> Result<()> {
    if let Some(name) = record.name() {
        write!(out, "{}\t", name).chain_err(|| "failed to write output row")?;
    }

    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
        String::from_utf8_lossy(&sequence),
        record.chrom(),
        record.start() as isize + site.start - offset + 1,
        record.start() as isize + site.end - offset,
        record.start() as isize + site.cutsite - offset + 1,
        site.strand.strand_symbol(),
        site.score,
    )
    .chain_err(|| "failed to write output row")
}

fn print_targets(
    index: &KMerIndex,
    idx: usize,
    record: &bed::Record,
    sequence: &[u8],
    out: &mut Box<dyn std::io::Write>,
) -> Result<()> {
    let enzyme = index.enzyme();
    let grna_len = enzyme.grna_len;
    let pam_len = enzyme.pam.len();
    let pam_offset = match enzyme.pam.position() {
        Position::Head => 0,
        Position::Tail => enzyme.grna_len - pam_len,
    };

    let num_sites = usize::max(sequence.len(), grna_len - 1) - grna_len - 1;
    let prefix = format!("{}. {}: ", idx + 1, record.name().unwrap_or(record.chrom()));
    let progress = progress::with_prefix(num_sites * 2, &prefix);

    let min_cutsite = u64::min(index.enzyme().grna_len as u64, record.start()) as isize;
    let max_cutsite = min_cutsite + (record.end() - record.start()) as isize;

    for site in collect_targets(index, sequence, &progress) {
        if site.cutsite >= min_cutsite && site.cutsite < max_cutsite {
            let mut target_seq = site.sequence.clone();
            for nuc in &mut target_seq[pam_offset..pam_offset + pam_len] {
                *nuc = nuc.to_ascii_lowercase();
            }

            print_target(record, &site, &target_seq, min_cutsite, out)?;
        }
    }

    progress.finish();

    Ok(())
}

fn collect_bed_targets(args: &FindArgs, index: &KMerIndex, bedfile: &str) -> Result<()> {
    // File handles are opened individually for better error reporting
    eprintln!("Finding target sites in {:?}", &args.targets);
    let fai = fasta::Index::with_fasta_file(&args.targets)
        .chain_err(|| format!("failed to open FASTA index file for {:?}", args.targets))?;
    let fasta = std::fs::File::open(&args.targets)
        .chain_err(|| format!("failed to open FASTA file {:?}", args.targets))?;
    let mut reader = fasta::IndexedReader::with_index(fasta, fai);

    let refseqs: HashMap<String, u64> = reader
        .index
        .sequences()
        .into_iter()
        .map(|seq| (seq.name, seq.len))
        .collect();

    eprintln!("  using regions from {:?}", bedfile);
    let mut beds = bed::Reader::from_file(bedfile)
        .chain_err(|| format!("failed to open BED file {:?}", bedfile))?;

    let mut out = open_file_or_stdout(&args.output)?;
    writeln!(
        out,
        "Region\tSequence\tContig\tStart\tEnd\tCutsite\tStrand\tScore"
    )
    .chain_err(|| "failed to write output header")?;

    for (idx, record) in beds.records().enumerate() {
        let record = record.chain_err(|| "failed to read BED record")?;
        let refseq_len = match refseqs.get(record.chrom()) {
            Some(&len) => len,
            None => return Err("foo".into()),
        };

        // Padding needed to find all cut-sites overlapping the target region
        let padding = index.enzyme().grna_len as isize;
        let start = isize::max(0, record.start() as isize - padding) as u64;
        let end = u64::min(refseq_len, record.end() + padding as u64);

        reader
            .fetch(record.chrom(), start, end)
            .chain_err(|| format!("failed to fetch {:?}", record))?;

        let mut sequence = Vec::new();
        reader
            .read(&mut sequence)
            .chain_err(|| format!("failed to read {:?}", record))?;

        sequence.make_ascii_uppercase();
        print_targets(index, idx, &record, &sequence, &mut out)?;
    }

    Ok(())
}

fn collect_all_targets(args: &FindArgs, index: &KMerIndex) -> Result<()> {
    eprintln!("Finding target sites in {:?}", &args.targets);
    let reader = fasta::Reader::from_file(&args.targets)
        .chain_err(|| format!("failed to open FASTA file {:?}", args.targets))?;

    let mut out = open_file_or_stdout(&args.output)?;
    writeln!(out, "Sequence\tContig\tStart\tEnd\tCutsite\tStrand\tScore")
        .chain_err(|| "failed to write output header")?;

    for (idx, target) in reader.records().enumerate() {
        let target = target.chain_err(|| "failed to read sequence")?;
        target.check().map_err(|v| ErrorKind::Msg(v.into()))?;

        let sequence = target.seq().to_ascii_uppercase();
        let mut record = bed::Record::new();
        record.set_chrom(target.id());
        record.set_end(sequence.len() as u64);

        print_targets(index, idx, &record, &sequence, &mut out)?;
    }

    Ok(())
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

    if let Some(bedfile) = &args.bedfile {
        collect_bed_targets(args, &index, bedfile)
    } else {
        collect_all_targets(args, &index)
    }
}
