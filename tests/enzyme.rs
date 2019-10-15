extern crate crispyr;
use crispyr::enzyme::Enzyme;

#[test]
fn test_enzyme_cas9() {
    assert_eq!(Enzyme::get("cas9").map(|v| v.name), Some("Cas9"));
    assert_eq!(Enzyme::get("Cas9").map(|v| v.name), Some("Cas9"));
    assert_eq!(Enzyme::get("cAs9").map(|v| v.name), Some("Cas9"));
    assert_eq!(Enzyme::get("CAS9").map(|v| v.name), Some("Cas9"));
}

#[test]
fn test_enzyme_mad7() {
    assert_eq!(Enzyme::get("mad7").map(|v| v.name), Some("Mad7"));
    assert_eq!(Enzyme::get("Mad7").map(|v| v.name), Some("Mad7"));
    assert_eq!(Enzyme::get("mAD7").map(|v| v.name), Some("Mad7"));
    assert_eq!(Enzyme::get("MAD7").map(|v| v.name), Some("Mad7"));
}

#[test]
fn test_enzyme_unknown() {
    assert_eq!(Enzyme::get("mad9"), None);
    assert_eq!(Enzyme::get("Cas7"), None);
    assert_eq!(Enzyme::get("Cat9"), None);
    assert_eq!(Enzyme::get("Foo"), None);
}
