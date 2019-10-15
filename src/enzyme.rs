use crate::pam::PAM;

#[derive(Clone, Debug, PartialEq)]
pub struct Enzyme {
    pub name: &'static str,
    pub extension: &'static str,

    pub grna_len: usize,

    pub pam: PAM,
    pub cutsite: Option<isize>,
}

impl Enzyme {
    pub fn get(name: &str) -> Option<Enzyme> {
        match name.to_ascii_lowercase().as_ref() {
            "cas9" => Some(Self::cas9()),
            "mad7" => Some(Self::mad7()),
            _ => None,
        }
    }

    pub fn cas9() -> Enzyme {
        Enzyme {
            name: "Cas9",
            extension: ".crispyr_cas9",

            grna_len: 23,

            pam: PAM::tail(b"NGG"),
            cutsite: Some(-3),
        }
    }

    pub fn mad7() -> Enzyme {
        Enzyme {
            name: "Mad7",
            extension: ".crispyr_mad7",

            grna_len: 25,

            pam: PAM::head(b"YTTN"),
            cutsite: None,
        }
    }
}
