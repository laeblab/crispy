use std::fmt::Debug;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::path::Path;

use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};

use crate::common::KMer;
use crate::constants::*;
use crate::enzyme::Enzyme;
use crate::errors::*;
use crate::pam;
use crate::progress;

const FORWARD_STRAND: u64 = 1 << 63;
const REVERSE_STRAND: u64 = 0;

#[derive(Debug, PartialEq, Clone)]
pub struct Position {
    refseq: u32,
    pos: i32,
    forward: bool,
}

impl Position {
    pub fn from_u64(raw: u64) -> Position {
        let forward = raw & FORWARD_STRAND != 0;
        let refseq = ((raw >> 32) & 0x7FFF_FFFF) as u32;
        let pos = (raw & 0xFFFF_FFFF) as i32;

        Position {
            refseq,
            pos,
            forward,
        }
    }

    pub fn to_u64(&self) -> u64 {
        assert!(
            self.refseq & 0x8000_0000 == 0,
            "cannot represent refseq using 31 bits"
        );

        let strand = if self.forward {
            FORWARD_STRAND
        } else {
            REVERSE_STRAND
        };
        let refseq = u64::from(self.refseq) << 32;
        let pos = u64::from(self.pos as u32);

        strand | refseq | pos
    }

    pub fn forward(refseq: u32, pos: i32) -> Position {
        Position {
            refseq,
            pos,
            forward: true,
        }
    }

    pub fn reverse(refseq: u32, pos: i32) -> Position {
        Position {
            refseq,
            pos,
            forward: false,
        }
    }

    pub fn refseq(&self) -> u32 {
        self.refseq
    }

    pub fn pos(&self) -> i32 {
        self.pos
    }

    pub fn strand(&self) -> char {
        if self.forward {
            '+'
        } else {
            '-'
        }
    }
}

pub enum KMerMap {
    Counts(Vec<u32>),
    Positions(Vec<Vec<Position>>),
}

impl KMerMap {
    pub fn len(&self) -> usize {
        match self {
            KMerMap::Counts(map) => map.iter().filter(|v| **v > 0).count(),
            KMerMap::Positions(map) => map.iter().filter(|v| !v.is_empty()).count(),
        }
    }

    pub fn add(&mut self, position: Position, kmer: KMer) {
        match self {
            KMerMap::Counts(map) => {
                map[kmer.key()] += 1;
            }
            KMerMap::Positions(map) => {
                map[kmer.key()].push(position);
            }
        }
    }
}

pub struct KMerIndex {
    enzyme: Enzyme,
    refseqs: Vec<String>,
    kmers: KMerMap,
}

impl KMerIndex {
    pub fn new(enzyme: &Enzyme, refseqs: Vec<String>, kmers: KMerMap) -> KMerIndex {
        KMerIndex {
            enzyme: enzyme.clone(),
            refseqs,
            kmers,
        }
    }

    pub fn refseqs(&self) -> &[String] {
        &self.refseqs
    }

    pub fn enzyme(&self) -> &Enzyme {
        &self.enzyme
    }

    pub fn kmer_count(&self) -> usize {
        self.kmers.len()
    }

    pub fn get_count(&self, kmer: KMer) -> Option<u32> {
        let count = match &self.kmers {
            KMerMap::Counts(map) => map[kmer.key()],
            KMerMap::Positions(map) => map[kmer.key()].len() as u32,
        };

        if count > 0 {
            Some(count)
        } else {
            None
        }
    }

    pub fn has_positions(&self) -> bool {
        match self.kmers {
            KMerMap::Counts(_) => false,
            KMerMap::Positions(_) => true,
        }
    }

    pub fn get_positions(&self, kmer: KMer) -> Option<&[Position]> {
        match &self.kmers {
            KMerMap::Counts(_) => None,
            KMerMap::Positions(map) => {
                let positions = &map[kmer.key()];

                if positions.is_empty() {
                    None
                } else {
                    Some(positions)
                }
            }
        }
    }

    pub fn read<P: AsRef<Path> + Debug>(filename: &P) -> Result<KMerIndex> {
        let file = File::open(filename).chain_err(|| "failed to open index file")?;
        let mut reader = BufReader::new(file);
        let mut buffer = Vec::new();

        reader
            .by_ref()
            .take(INDEX_HEADER.len() as u64)
            .read_to_end(&mut buffer)
            .chain_err(|| "failed to read index header")?;
        if buffer != INDEX_HEADER {
            println!("{:?}", buffer);
            return Err("file is not a valid CRISPyR index file".into());
        }

        let version = reader
            .read_u8()
            .chain_err(|| "failed to read index version number")?;
        if version < INDEX_VERSION {
            return Err("index file is outdated; please re-index genome".into());
        } else if version > INDEX_VERSION {
            return Err("index generated using newer version of CRISPyR;
                        please upgrade CRISPyR or re-index genome"
                .into());
        }

        buffer.clear();
        let enzyme_len = reader
            .read_u8()
            .chain_err(|| "failed to read length of enzyme name")?;
        reader
            .by_ref()
            .take(u64::from(enzyme_len))
            .read_to_end(&mut buffer)
            .chain_err(|| "failed to read enzyme name")?;

        let enzyme_name =
            std::str::from_utf8(&buffer).chain_err(|| "failed to decode enzyme name")?;

        let enzyme = match Enzyme::get(&enzyme_name) {
            Some(value) => value,
            None => return Err(format!("unknown enzyme {:?}", enzyme_name).into()),
        };

        let flags = reader
            .read_u64::<LittleEndian>()
            .chain_err(|| "failed to read index flags")?;

        let (refseqs, kmers) = if flags & FLAG_POSITIONS != 0 {
            (
                Self::read_refseqs(&mut reader)?,
                Self::read_positions(&mut reader)?,
            )
        } else {
            (Vec::new(), Self::read_counts(&mut reader)?)
        };

        Ok(KMerIndex {
            enzyme,
            refseqs,
            kmers,
        })
    }

    pub fn write<P: AsRef<Path> + Debug>(&self, filename: P) -> Result<()> {
        let file = File::create(filename).chain_err(|| "failed to create index file")?;
        let mut writer = BufWriter::new(file);

        writer
            .write(INDEX_HEADER)
            .chain_err(|| "failed to write index header")?;
        writer
            .write_u8(INDEX_VERSION)
            .chain_err(|| "failed to write index version")?;

        let name = self.enzyme.name.as_bytes();
        writer
            .write_u8(name.len() as u8)
            .chain_err(|| "failed to write enzyme name length")?;
        writer
            .write(&name)
            .chain_err(|| "failed to write enzyme name")?;

        match &self.kmers {
            KMerMap::Counts(map) => {
                writer
                    .write_u64::<LittleEndian>(0)
                    .chain_err(|| "failed to write index flags")?;

                writer
                    .write_u64::<LittleEndian>(self.kmers.len() as u64)
                    .chain_err(|| "failed to write number of unique kmers")?;
                Self::write_counts(&mut writer, map).chain_err(|| "failed to write KMers")
            }
            KMerMap::Positions(map) => {
                writer
                    .write_u64::<LittleEndian>(FLAG_POSITIONS)
                    .chain_err(|| "failed to write index flags")?;

                Self::write_refseqs(&mut writer, &self.refseqs)
                    .chain_err(|| "failed to write reference sequenec names")?;

                writer
                    .write_u64::<LittleEndian>(self.kmers.len() as u64)
                    .chain_err(|| "failed to write number of unique kmers")?;
                Self::write_positions(&mut writer, map).chain_err(|| "failed to write KMers")
            }
        }
    }

    fn read_counts(reader: &mut BufReader<File>) -> Result<KMerMap> {
        let kmer_count = reader
            .read_u64::<LittleEndian>()
            .chain_err(|| "failed to read number of unique sites")?;
        let mut kmers = vec![0; KMER_COUNT];

        let progress = progress::default(kmer_count as usize);
        for _ in 0..kmer_count {
            let kmer = reader
                .read_u32::<LittleEndian>()
                .chain_err(|| "failed to read K-mer")?;
            let count = reader
                .read_u32::<LittleEndian>()
                .chain_err(|| "failed to read K-mer frequency")?;

            kmers[kmer as usize] += count;
            progress.inc(1);
        }

        progress.finish();

        Ok(KMerMap::Counts(kmers))
    }

    fn write_counts(writer: &mut BufWriter<File>, map: &[u32]) -> Result<()> {
        for (kmer, count) in map.iter().enumerate() {
            if *count > 0 {
                writer
                    .write_u32::<LittleEndian>(kmer as u32)
                    .chain_err(|| "failed to write kmer")?;
                writer
                    .write_u32::<LittleEndian>(*count)
                    .chain_err(|| "failed to write kmer count")?;
            }
        }

        Ok(())
    }

    fn write_refseqs(writer: &mut BufWriter<File>, refseqs: &[String]) -> Result<()> {
        writer
            .write_u64::<LittleEndian>(refseqs.len() as u64)
            .chain_err(|| "failed to write number of reference sequences")?;

        for name in refseqs {
            let bytes = name.as_bytes();

            writer
                .write_u16::<LittleEndian>(bytes.len() as u16)
                .chain_err(|| "failed to write length of reference sequence name")?;
            writer
                .write(bytes)
                .chain_err(|| "failed to write reference sequence name")?;
        }

        Ok(())
    }

    fn write_positions(writer: &mut BufWriter<File>, map: &[Vec<Position>]) -> Result<()> {
        for (kmer, positions) in map.iter().enumerate() {
            if !positions.is_empty() {
                writer
                    .write_u32::<LittleEndian>(kmer as u32)
                    .chain_err(|| "failed to write kmer")?;
                writer
                    .write_u32::<LittleEndian>(positions.len() as u32)
                    .chain_err(|| "failed to write number of hits for kmer")?;

                for position in positions {
                    writer
                        .write_u64::<LittleEndian>(position.to_u64())
                        .chain_err(|| "failed to write kmer position")?;
                }
            }
        }

        Ok(())
    }

    fn read_refseqs(reader: &mut BufReader<File>) -> Result<Vec<String>> {
        let refseq_count = reader
            .read_u64::<LittleEndian>()
            .chain_err(|| "failed to read number of reference sequences")?;

        let progress = progress::default(refseq_count as usize);
        let mut refseqs = Vec::with_capacity(refseq_count as usize);
        for _ in 0..refseq_count {
            let refseq_len = reader
                .read_u16::<LittleEndian>()
                .chain_err(|| "failed to read reference sequence length")?;

            let mut buffer = vec![0; refseq_len as usize];
            reader
                .read_exact(&mut buffer)
                .chain_err(|| "failed to read reference sequence name")?;

            let name = String::from_utf8(buffer)
                .chain_err(|| "failed to parse reference sequence name")?;

            refseqs.push(name);
            progress.inc(1);
        }

        progress.finish();

        Ok(refseqs)
    }

    fn read_positions(reader: &mut BufReader<File>) -> Result<KMerMap> {
        let mut kmers = vec![Vec::new(); KMER_COUNT];
        let kmer_count = reader
            .read_u64::<LittleEndian>()
            .chain_err(|| "failed to read number of reference sequences")?;

        let progress = progress::default(kmer_count as usize);
        for _ in 0..kmer_count {
            let kmer = reader
                .read_u32::<LittleEndian>()
                .chain_err(|| "failed to read K-mer")?;
            let count = reader
                .read_u32::<LittleEndian>()
                .chain_err(|| "failed to read K-mer frequency")?;

            let mut positions = Vec::new();
            for _ in 0..count {
                let refseq = reader
                    .read_u64::<LittleEndian>()
                    .chain_err(|| "failed to read kmer position")?;

                positions.push(Position::from_u64(refseq));
            }

            kmers[kmer as usize] = positions;
            progress.inc(1);
        }

        progress.finish();

        Ok(KMerMap::Positions(kmers))
    }

    pub fn summarize(&self) -> String {
        let end = match self.enzyme.pam.position() {
            pam::Position::Head => "5'",
            pam::Position::Tail => "3'",
        };

        format!(
            "Index contains {} unique K-mers for {} with {} PAM sequence {}",
            self.kmer_count(),
            self.enzyme.name,
            end,
            self.enzyme.pam.to_string()
        )
    }
}
