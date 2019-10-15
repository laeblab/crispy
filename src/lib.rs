// error_chain macro can recurse deeply
#![recursion_limit = "1024"]

#[macro_use]
extern crate error_chain;
#[macro_use(lazy_static)]
extern crate lazy_static;

pub mod args;
pub mod commands;
pub mod common;
pub mod constants;
pub mod enzyme;
pub mod errors;
pub mod index;
pub mod iupac;
pub mod pam;
pub mod progress;
pub mod score;
pub mod table;
