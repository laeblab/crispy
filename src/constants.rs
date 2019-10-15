pub const KMER_LEN: usize = 13;

pub const KMER_COUNT: usize = 2 << (2 * KMER_LEN - 1);

pub const INDEX_HEADER: &[u8] = b"CRISPyR";
pub const INDEX_VERSION: u8 = 4;

// Index flag indicating if PAM positions have been saved
pub const FLAG_POSITIONS: u64 = 0b1;
