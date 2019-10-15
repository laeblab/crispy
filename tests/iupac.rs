extern crate crispyr;

use crispyr::iupac::matches;

#[test]
fn test_basic_nucleotides() {
    // first
    assert!(matches(b'A', b'A'));
    assert!(!matches(b'A', b'C'));
    assert!(!matches(b'A', b'G'));
    assert!(!matches(b'A', b'T'));
    // last
    assert!(!matches(b'T', b'A'));
    assert!(!matches(b'T', b'C'));
    assert!(!matches(b'T', b'G'));
    assert!(matches(b'T', b'T'));
}

#[test]
fn test_degenerate_nucleotide() {
    assert!(!matches(b'N', b'@'));
    assert!(matches(b'N', b'A')); // A
    assert!(matches(b'N', b'B')); // C, G, T
    assert!(matches(b'N', b'C')); // C
    assert!(matches(b'N', b'D')); // A, G, T
    assert!(!matches(b'N', b'E'));
    assert!(!matches(b'N', b'['));
}

#[test]
fn test_non_nucleotide_values() {
    assert!(matches(b'I', b'I'));
    assert!(!matches(b'I', b'J'));
    assert!(matches(b'!', b'!'));
    assert!(!matches(b'!', b'?'));
}
