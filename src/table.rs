use std::fmt::Debug;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

use crate::errors::*;

pub fn read<P: AsRef<Path> + Debug>(path: &P) -> Result<Vec<Vec<String>>> {
    let file = File::open(path).chain_err(|| "failed to open table")?;
    let reader = BufReader::new(file);

    let mut table: Vec<Vec<String>> = Vec::new();
    for line in reader.lines() {
        let line = line.chain_err(|| "error reading line from table")?;
        if !line.trim().is_empty() {
            table.push(line.split('\t').map(|v| v.to_string()).collect());
        }
    }

    Ok(table)
}
