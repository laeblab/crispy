// error_chain macro can recurse deeply
#![recursion_limit = "1024"]

#[macro_use]
extern crate error_chain;
#[macro_use(lazy_static)]
extern crate lazy_static;

mod args;
mod commands;
mod common;
mod constants;
mod enzyme;
mod errors;
mod index;
mod iupac;
mod pam;
mod progress;
mod score;
mod table;

fn print_err(e: &errors::Error) {
    use error_chain::ChainedError;
    use std::io::Write; // trait which holds `display_chain`
    let stderr = &mut ::std::io::stderr();
    let errmsg = "Error writing to stderr";

    writeln!(stderr, "{}", e.display_chain()).expect(errmsg);
}

fn inner_main() -> errors::Result<()> {
    match args::parse_args()? {
        args::Args::Find(args) => commands::find::main(&args),
        args::Args::Index(args) => commands::index::main(&args),
        args::Args::OffTargets(args) => commands::offtargets::main(&args),
        args::Args::Score(args) => commands::score::main(&args),
        args::Args::None => Ok(()),
    }
}

fn main() {
    if let Err(e) = inner_main() {
        print_err(&e);

        ::std::process::exit(1);
    } else {
        ::std::process::exit(0);
    }
}
