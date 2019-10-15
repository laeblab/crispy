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
    pub threads: usize,
}

#[derive(Debug)]
pub struct FindArgs {
    pub index: String,
    pub targets: String,
    pub threads: usize,
}

#[derive(Debug)]
pub struct OffTargetsArgs {
    pub index: String,
    pub table: String,
}

pub enum Args {
    Index(IndexArgs),
    Score(ScoreArgs),
    Find(FindArgs),
    OffTargets(OffTargetsArgs),
    None,
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
        .arg(
            Arg::with_name("index")
                .help("Path to CRISPyR index file.")
                .required(true),
        )
        .arg(
            Arg::with_name("table")
                .help("Table containing target sequences.")
                .required(true),
        )
        .arg(
            Arg::with_name("threads")
                .long("threads")
                .takes_value(true)
                .allow_hyphen_values(true)
                .number_of_values(1)
                .default_value("0")
                .help("Number of threads used for computation (0 for automatic)."),
        )
}

fn find_command<'a, 'b>() -> App<'a, 'b> {
    SubCommand::with_name("find")
        .about("Find and score gRNA targets in FASTA sequence(s)")
        .arg(
            Arg::with_name("index")
                .help("Path to CRISPyR index file.")
                .required(true),
        )
        .arg(
            Arg::with_name("targets")
                .help("FASTA file containing one or more sequences.")
                .required(true),
        )
        .arg(
            Arg::with_name("threads")
                .long("threads")
                .takes_value(true)
                .allow_hyphen_values(true)
                .number_of_values(1)
                .default_value("0")
                .help("Number of threads used for computation (0 for automatic)."),
        )
}

fn off_targets_command<'a, 'b>() -> App<'a, 'b> {
    SubCommand::with_name("offtargets")
        .about("Print table of off targets for each target sequence")
        .arg(
            Arg::with_name("index")
                .help("Path to CRISPyR index file.")
                .required(true),
        )
        .arg(
            Arg::with_name("table")
                .help("Table containing target sequences.")
                .required(true),
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

pub fn parse_args() -> Result<Args> {
    let matches = App::new("CRISPyR")
        .version("0.0.1")
        .author("Mikkel Schubert")
        .subcommand(index_command())
        .subcommand(score_command())
        .subcommand(find_command())
        .subcommand(off_targets_command())
        .get_matches();

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
            threads: parse_threads(matches)?,
        }))
    } else if let Some(matches) = matches.subcommand_matches("offtargets") {
        Ok(Args::OffTargets(OffTargetsArgs {
            index: get_string(matches, "index")?,
            table: get_string(matches, "table")?,
        }))
    } else if let Some(matches) = matches.subcommand_matches("find") {
        Ok(Args::Find(FindArgs {
            index: get_string(matches, "index")?,
            targets: get_string(matches, "targets")?,
            threads: parse_threads(matches)?,
        }))
    } else {
        eprintln!("{}", matches.usage());

        Ok(Args::None)
    }
}
