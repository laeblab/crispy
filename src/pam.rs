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
            match self.position {
                Position::Head => {
                    encode_dna(&window[self.len()..self.len() + KMER_LEN]).map(|v| (0, v))
                }
                Position::Tail => {
                    let pam_start = window.len() - self.len();
                    let kmer_start = pam_start - KMER_LEN;
                    let kmer = encode_dna(&window[kmer_start..pam_start]);

                    kmer.map(|v| (pam_start, v))
                }
            }
        } else {
            None
        }
    }

    pub fn len(&self) -> usize {
        self.pam.len()
    }

    pub fn to_string(&self) -> Cow<str> {
        String::from_utf8_lossy(&self.pam)
    }
}
