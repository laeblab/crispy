use crate::constants::KMER_LEN;

#[derive(Hash, PartialEq, Eq, Debug, Clone, Copy)]
pub struct KMer(pub u32);

impl KMer {
    pub fn new(kmer: u32) -> KMer {
        KMer(kmer)
    }
}

impl KMer {
    pub fn key(self) -> usize {
        self.0 as usize
    }
}

pub fn encode_dna(seq: &[u8]) -> Option<KMer> {
    assert!(seq.len() == KMER_LEN);

    let mut encoded_dna = 0;
    for nuc in seq {
        encoded_dna = (encoded_dna << 2)
            | match *nuc {
                b'a' | b'A' => 0,
                b'c' | b'C' => 1,
                b'g' | b'G' => 2,
                b't' | b'T' => 3,
                _ => return None,
            };
    }

    Some(KMer::new(encoded_dna))
}
