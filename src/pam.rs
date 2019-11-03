use std::borrow::Cow;

use crate::common::{encode_dna, KMer};
use crate::constants::KMER_LEN;
use crate::iupac;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Position {
    Head,
    Tail,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PAM {
    pam: Vec<u8>,
    position: Position,
}

impl PAM {
    pub fn head(query: &[u8]) -> PAM {
        PAM {
            pam: query.to_owned(),
            position: Position::Head,
        }
    }

    pub fn tail(query: &[u8]) -> PAM {
        PAM {
            pam: query.to_owned(),
            position: Position::Tail,
        }
    }

    pub fn position(&self) -> Position {
        self.position
    }

    pub fn matches(&self, window: &[u8]) -> bool {
        if self.len() <= window.len() {
            let iupac_wrapper = |(&query, &candidate)| iupac::matches(query, candidate);

            match self.position {
                Position::Head => self.pam.iter().zip(window.iter()).all(iupac_wrapper),
                Position::Tail => self
                    .pam
                    .iter()
                    .rev()
                    .zip(window.iter().rev())
                    .all(iupac_wrapper),
            }
        } else {
            false
        }
    }

    pub fn kmer(&self, window: &[u8]) -> Option<(usize, KMer)> {
        if self.pam.len() + KMER_LEN <= window.len() && self.matches(window) {
            let kmer = encode_dna(self.kmer_slice(window));

            match self.position {
                Position::Head => kmer.map(|v| (0, v)),
                Position::Tail => kmer.map(|v| (window.len() - self.len(), v)),
            }
        } else {
            None
        }
    }

    /// Returns slice containing kmer. The string may be any length, but is assumed to start with
    /// the PAM for 5' PAMs and to end with the PAM for 3' PAMs.
    pub fn kmer_slice<'a>(&self, seq: &'a [u8]) -> &'a [u8] {
        let (start, end) = self.kmer_pos(seq);

        &seq[start as usize..end as usize]
    }

    /// Returns mutable slice containing kmer. The string may be any length, but is assumed to
    /// start with the PAM for 5' PAMs and to end with the PAM for 3' PAMs.
    pub fn kmer_slice_mut<'a>(&self, seq: &'a mut [u8]) -> &'a mut [u8] {
        let (start, end) = self.kmer_pos(seq);

        &mut seq[start..end]
    }

    fn kmer_pos(&self, seq: &[u8]) -> (usize, usize) {
        let seq_len = seq.len() as isize;
        let pam_len = self.pam.len() as isize;
        let kmer_len = KMER_LEN as isize;

        match self.position {
            Position::Head => (
                isize::min(seq_len, pam_len) as usize,
                isize::min(seq_len, pam_len + kmer_len) as usize,
            ),
            Position::Tail => (
                isize::max(0, seq_len - pam_len - kmer_len) as usize,
                isize::max(0, seq_len - pam_len) as usize,
            ),
        }
    }

    pub fn len(&self) -> usize {
        self.pam.len()
    }

    pub fn to_string(&self) -> Cow<str> {
        String::from_utf8_lossy(&self.pam)
    }
}
