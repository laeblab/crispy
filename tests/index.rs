extern crate crispyr;

use crispyr::index::*;

#[test]
fn test_position_to_u64_forward_zeroes() {
    let pos = Position::forward(0, 0);
    let repr = pos.to_u64();

    assert_eq!(pos, Position::from_u64(repr));
}

#[test]
fn test_position_to_u64_reverse_zeroes() {
    let pos = Position::reverse(0, 0);
    let repr = pos.to_u64();

    assert_eq!(pos, Position::from_u64(repr));
}

#[test]
fn test_position_to_u64_forward_ref_no_pos() {
    let pos = Position::forward(274, 0);
    let repr = pos.to_u64();

    assert_eq!(pos, Position::from_u64(repr));
}

#[test]
fn test_position_to_u64_reverse_ref_no_pos() {
    let pos = Position::reverse(274, 0);
    let repr = pos.to_u64();

    assert_eq!(pos, Position::from_u64(repr));
}

#[test]
fn test_position_to_u64_forward_pos_no_ref() {
    let pos = Position::forward(0, 7913);
    let repr = pos.to_u64();

    assert_eq!(pos, Position::from_u64(repr));
}

#[test]
fn test_position_to_u64_reverse_pos_no_ref() {
    let pos = Position::reverse(0, 7913);
    let repr = pos.to_u64();

    assert_eq!(pos, Position::from_u64(repr));
}

#[test]
fn test_position_to_u64_forward_neg_pos_no_ref() {
    let pos = Position::forward(0, -7913);
    let repr = pos.to_u64();

    assert_eq!(pos, Position::from_u64(repr));
}

#[test]
fn test_position_to_u64_reverse_neg_pos_no_ref() {
    let pos = Position::reverse(0, -7913);
    let repr = pos.to_u64();

    assert_eq!(pos, Position::from_u64(repr));
}

#[test]
fn test_position_to_u64_forward_pos_and_ref() {
    let pos = Position::forward(17, 7913);
    let repr = pos.to_u64();

    assert_eq!(pos, Position::from_u64(repr));
}

#[test]
fn test_position_to_u64_reverse_pos_and_ref() {
    let pos = Position::reverse(17, 7913);
    let repr = pos.to_u64();

    assert_eq!(pos, Position::from_u64(repr));
}

#[test]
fn test_position_to_u64_forward_neg_pos_and_ref() {
    let pos = Position::forward(17, -7913);
    let repr = pos.to_u64();

    assert_eq!(pos, Position::from_u64(repr));
}

#[test]
fn test_position_to_u64_reverse_neg_pos_and_ref() {
    let pos = Position::reverse(17, -7913);
    let repr = pos.to_u64();

    assert_eq!(pos, Position::from_u64(repr));
}
