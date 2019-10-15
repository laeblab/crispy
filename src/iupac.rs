const MAPPING: [(u8, &[u8]); 15] = [
    (b'A', b"A"),
    (b'C', b"C"),
    (b'G', b"G"),
    (b'T', b"T"),
    (b'R', b"AGR"),
    (b'Y', b"CTY"),
    (b'S', b"GCS"),
    (b'W', b"ATW"),
    (b'K', b"GTK"),
    (b'M', b"ACM"),
    (b'B', b"CGTB"),
    (b'D', b"AGTD"),
    (b'H', b"ACTH"),
    (b'V', b"ACGV"),
    (b'N', b"ACGTRYSWKMBDHVN"),
];

lazy_static! {
    static ref IUPAC: Vec<bool> = {
        let mut table = vec![false; 26 * 26];

        for (query, matches) in &MAPPING {
            for &candidate in matches.iter() {
                table[offset(*query, candidate)] = true;
            }
        }

        table
    };
}

fn offset(query: u8, candidate: u8) -> usize {
    (query - b'A') as usize * 26 + (candidate - b'A') as usize
}

pub fn matches(query: u8, candidate: u8) -> bool {
    if query == candidate {
        true
    } else if b'A' <= query && query <= b'Z' && b'A' <= candidate && candidate <= b'Z' {
        IUPAC[offset(query, candidate)]
    } else {
        false
    }
}
