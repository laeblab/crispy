use std::io;
use std::io::Write;

use clap::{App, Arg, ArgMatches, SubCommand};

use crate::enzyme::Enzyme;
use crate::errors::*;

#[derive(Debug)]
pub struct IndexArgs {
    pub fasta: String,
    pub output: Option<String>,
    pub enzyme: Enzyme,
    pub positions: bool,
}

#[derive(Debug)]
pub struct ScoreArgs {
    pub index: String,
    pub table: String,
    pub output: Option<String>,
    pub threads: usize,
}

#[derive(Debug)]
pub struct FindArgs {
    pub index: String,
    pub targets: String,
    pub output: Option<String>,
    pub bedfile: Option<String>,
    pub threads: usize,
}

#[derive(Debug)]
pub struct OffTargetsArgs {
    pub index: String,
    pub table: String,
    pub fasta: Option<String>,
    pub output: Option<String>,
    pub min_score: u64,
}

pub enum Args {
    Index(IndexArgs),
    Score(ScoreArgs),
    Find(FindArgs),
    OffTargets(OffTargetsArgs),
    None,
}

/// Command-line parameter for CRISPyR index
fn args_index<'a, 'b>() -> Arg<'a, 'b> {
    Arg::with_name("index")
        .help("Path to CRISPyR index file.")
        .required(true)
}

/// Command-line option for specifying output files
fn args_output<'a, 'b>() -> Arg<'a, 'b> {
    Arg::with_name("output")
        .long("output")
        .short("o")
        .takes_value(true)
        .number_of_values(1)
        .help("Write output to file instead of STDOUT.")
}

/// Command-line option for specifying the number of threads used
fn args_threads<'a, 'b>() -> Arg<'a, 'b> {
    Arg::with_name("threads")
        .long("threads")
        .takes_value(true)
        .allow_hyphen_values(true)
        .number_of_values(1)
        .default_value("0")
        .help("Number of threads used for computation (0 for automatic).")
}

fn index_command<'a, 'b>() -> App<'a, 'b> {
    SubCommand::with_name("index")
        .about("Index genome in FASTA format")
        .arg(Arg::with_name("fasta").help("TODO").required(true))
        .arg(Arg::with_name("output"))
        .arg(
            Arg::with_name("enzyme")
                .long("enzyme")
                .takes_value(true)
                .default_value("Cas9")
                .help("Endonuclease enzyme used; either Cas9 or Mad7."),
        )
        .arg(
            Arg::with_name("positions")
                .long("positions")
                .help("Save cut-site positions (greatly increases index size)"),
        )
}

fn score_command<'a, 'b>() -> App<'a, 'b> {
    SubCommand::with_name("score")
        .about("Score table of gRNA targets using indexed genome")
        .arg(args_index())
        .arg(
            Arg::with_name("table")
                .help("Table containing target sequences.")
                .required(true),
        )
        .arg(args_output())
        .arg(args_threads())
}

fn find_command<'a, 'b>() -> App<'a, 'b> {
    SubCommand::with_name("find")
        .about("Find and score gRNA targets in FASTA sequence(s)")
        .arg(args_index())
        .arg(
            Arg::with_name("targets")
                .help("FASTA file containing one or more sequences.")
                .required(true),
        )
        .arg(
            Arg::with_name("bedfile")
                .long("bed")
                .takes_value(true)
                .number_of_values(1)
                .help(
                    "Collect gRNAs with cut-sites inside the named regions in this \
                     BED file. Note that cut-sites may not be exact and that it is \
                     therefore safest to pick gRNAs with cut-sites some number of bp \
                     from the start and the end of the target region. Also note that \
                     overlapping regions are not merged.",
                ),
        )
        .arg(args_output())
        .arg(args_threads())
}

fn off_targets_command<'a, 'b>() -> App<'a, 'b> {
    SubCommand::with_name("offtargets")
        .about("Print table of off targets for each target sequence")
        .arg(args_index())
        .arg(
            Arg::with_name("table")
                .help("Table containing target sequences.")
                .required(true),
        )
        .arg(
            Arg::with_name("fasta")
                .help(
                    "FASTA file for index; if not specified, CRISPyR will \
                     try using the index filename without the extension.",
                )
                .required(false),
        )
        .arg(args_output())
        .arg(
            Arg::with_name("min_score")
                .long("min-score")
                .takes_value(true)
                .number_of_values(1)
                .default_value("0")
                .help("Minimum score of off-target (500 for exact kmer matches)"),
        )
        .alias("off_targets")
}

fn get_str<'a>(matches: &'a ArgMatches, key: &str) -> Result<&'a str> {
    match matches.value_of(key) {
        Some(value) => Ok(value),
        None => Err(format!("Required option {:?} not set", key).into()),
    }
}

fn get_string(matches: &ArgMatches, key: &str) -> Result<String> {
    get_str(matches, key).map(|v| v.into())
}

fn parse_threads(matches: &ArgMatches) -> Result<usize> {
    let s = get_str(matches, "threads")?;

    match usize::from_str_radix(s, 10) {
        Ok(v) => Ok(v),
        Err(err) => Err(format!("Invalid --threads ({:?}) value: {}", s, err).into()),
    }
}

fn parse_min_score(matches: &ArgMatches) -> Result<u64> {
    let s = get_str(matches, "min_score")?;

    match u64::from_str_radix(s, 10) {
        Ok(v) => Ok(v),
        Err(err) => Err(format!("Invalid --min-score ({:?}) value: {}", s, err).into()),
    }
}

fn new_parser<'a, 'b>() -> App<'a, 'b> {
    App::new("CRISPyR")
        .version("0.2.0")
        .author("Mikkel Schubert")
        .subcommand(index_command())
        .subcommand(score_command())
        .subcommand(find_command())
        .subcommand(off_targets_command())
}

pub fn parse_args() -> Result<Args> {
    let matches = new_parser().get_matches();

    if let Some(matches) = matches.subcommand_matches("index") {
        let enzyme_str = get_str(matches, "enzyme")?;
        let enzyme = match Enzyme::get(enzyme_str) {
            Some(enzyme) => enzyme,
            None => return Err(format!("Unknown enzyme {:?}", enzyme_str).into()),
        };

        Ok(Args::Index(IndexArgs {
            fasta: get_string(matches, "fasta")?,
            output: matches.value_of("output").map(|s| s.to_string()),
            enzyme,
            positions: matches.is_present("positions"),
        }))
    } else if let Some(matches) = matches.subcommand_matches("score") {
        Ok(Args::Score(ScoreArgs {
            index: get_string(matches, "index")?,
            table: get_string(matches, "table")?,
            output: matches.value_of("output").map(|s| s.to_string()),
            threads: parse_threads(matches)?,
        }))
    } else if let Some(matches) = matches.subcommand_matches("offtargets") {
        Ok(Args::OffTargets(OffTargetsArgs {
            index: get_string(matches, "index")?,
            table: get_string(matches, "table")?,
            fasta: matches.value_of("fasta").map(|s| s.to_string()),
            output: matches.value_of("output").map(|s| s.to_string()),
            min_score: parse_min_score(matches)?,
        }))
    } else if let Some(matches) = matches.subcommand_matches("find") {
        Ok(Args::Find(FindArgs {
            index: get_string(matches, "index")?,
            targets: get_string(matches, "targets")?,
            bedfile: matches.value_of("bedfile").map(|s| s.to_string()),
            output: matches.value_of("output").map(|s| s.to_string()),
            threads: parse_threads(matches)?,
        }))
    } else {
        let mut out = io::stderr();
        let _ = new_parser().write_help(&mut out);
        let _ = writeln!(out);

        Ok(Args::None)
    }
}
