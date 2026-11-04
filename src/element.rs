#[derive(Clone, Copy, Debug)]
pub struct Properties {
    pub symbol: &'static str,
    pub name: &'static str,
    pub weight: f64,
    pub covalent_radius: f64,
    /// Van der Waals radius in Ångström (Alvarez, 2013), for non-bonded
    /// contact distance. `None` for elements Alvarez's survey has no
    /// experimentally-derived value for: promethium, polonium through radium,
    /// fermium onward, and every element past lawrencium — all short-lived
    /// synthetic elements with no measured van der Waals data.
    pub van_der_waals_radius: Option<f64>,
    /// Allowed valences of the neutral atom: the sums of bond orders it may have (RDKit's
    /// table). Phosphorus is `[3, 5]` and sulfur `[2, 4, 6]` because each has more than one
    /// normal bonding state; noble gases are `[0]`. Empty means no constraint: transition
    /// metals, lanthanides, actinides and the alkali and alkaline-earth metals, whose bonding
    /// in structures is ionic or coordinative and not described by a fixed valence.
    pub valences: &'static [u8],
}

#[repr(u8)]
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub enum Element {
    H = 1,
    He,
    Li,
    Be,
    B,
    C,
    N,
    O,
    F,
    Ne,
    Na,
    Mg,
    Al,
    Si,
    P,
    S,
    Cl,
    Ar,
    K,
    Ca,
    Sc,
    Ti,
    V,
    Cr,
    Mn,
    Fe,
    Co,
    Ni,
    Cu,
    Zn,
    Ga,
    Ge,
    As,
    Se,
    Br,
    Kr,
    Rb,
    Sr,
    Y,
    Zr,
    Nb,
    Mo,
    Tc,
    Ru,
    Rh,
    Pd,
    Ag,
    Cd,
    In,
    Sn,
    Sb,
    Te,
    I,
    Xe,
    Cs,
    Ba,
    La,
    Ce,
    Pr,
    Nd,
    Pm,
    Sm,
    Eu,
    Gd,
    Tb,
    Dy,
    Ho,
    Er,
    Tm,
    Yb,
    Lu,
    Hf,
    Ta,
    W,
    Re,
    Os,
    Ir,
    Pt,
    Au,
    Hg,
    Tl,
    Pb,
    Bi,
    Po,
    At,
    Rn,
    Fr,
    Ra,
    Ac,
    Th,
    Pa,
    U,
    Np,
    Pu,
    Am,
    Cm,
    Bk,
    Cf,
    Es,
    Fm,
    Md,
    No,
    Lr,
    Rf,
    Db,
    Sg,
    Bh,
    Hs,
    Mt,
    Ds,
    Rg,
    Cn,
    Nh,
    Fl,
    Mc,
    Lv,
    Ts,
    Og,
}

impl Element {
    pub const fn atomic_number(self) -> u8 {
        self as u8
    }

    pub const fn symbol(self) -> &'static str {
        self.properties().symbol
    }

    pub const fn name(self) -> &'static str {
        self.properties().name
    }

    pub const fn weight(self) -> f64 {
        self.properties().weight
    }

    pub const fn covalent_radius(self) -> f64 {
        self.properties().covalent_radius
    }

    pub const fn van_der_waals_radius(self) -> Option<f64> {
        self.properties().van_der_waals_radius
    }

    /// Allowed valences of the neutral atom; see [`Properties::valences`].
    pub const fn valences(self) -> &'static [u8] {
        self.properties().valences
    }

    pub fn from_symbol(symbol: &str) -> Option<Element> {
        match symbol {
            "H" => Some(Element::H),
            "He" => Some(Element::He),
            "Li" => Some(Element::Li),
            "Be" => Some(Element::Be),
            "B" => Some(Element::B),
            "C" => Some(Element::C),
            "N" => Some(Element::N),
            "O" => Some(Element::O),
            "F" => Some(Element::F),
            "Ne" => Some(Element::Ne),
            "Na" => Some(Element::Na),
            "Mg" => Some(Element::Mg),
            "Al" => Some(Element::Al),
            "Si" => Some(Element::Si),
            "P" => Some(Element::P),
            "S" => Some(Element::S),
            "Cl" => Some(Element::Cl),
            "Ar" => Some(Element::Ar),
            "K" => Some(Element::K),
            "Ca" => Some(Element::Ca),
            "Sc" => Some(Element::Sc),
            "Ti" => Some(Element::Ti),
            "V" => Some(Element::V),
            "Cr" => Some(Element::Cr),
            "Mn" => Some(Element::Mn),
            "Fe" => Some(Element::Fe),
            "Co" => Some(Element::Co),
            "Ni" => Some(Element::Ni),
            "Cu" => Some(Element::Cu),
            "Zn" => Some(Element::Zn),
            "Ga" => Some(Element::Ga),
            "Ge" => Some(Element::Ge),
            "As" => Some(Element::As),
            "Se" => Some(Element::Se),
            "Br" => Some(Element::Br),
            "Kr" => Some(Element::Kr),
            "Rb" => Some(Element::Rb),
            "Sr" => Some(Element::Sr),
            "Y" => Some(Element::Y),
            "Zr" => Some(Element::Zr),
            "Nb" => Some(Element::Nb),
            "Mo" => Some(Element::Mo),
            "Tc" => Some(Element::Tc),
            "Ru" => Some(Element::Ru),
            "Rh" => Some(Element::Rh),
            "Pd" => Some(Element::Pd),
            "Ag" => Some(Element::Ag),
            "Cd" => Some(Element::Cd),
            "In" => Some(Element::In),
            "Sn" => Some(Element::Sn),
            "Sb" => Some(Element::Sb),
            "Te" => Some(Element::Te),
            "I" => Some(Element::I),
            "Xe" => Some(Element::Xe),
            "Cs" => Some(Element::Cs),
            "Ba" => Some(Element::Ba),
            "La" => Some(Element::La),
            "Ce" => Some(Element::Ce),
            "Pr" => Some(Element::Pr),
            "Nd" => Some(Element::Nd),
            "Pm" => Some(Element::Pm),
            "Sm" => Some(Element::Sm),
            "Eu" => Some(Element::Eu),
            "Gd" => Some(Element::Gd),
            "Tb" => Some(Element::Tb),
            "Dy" => Some(Element::Dy),
            "Ho" => Some(Element::Ho),
            "Er" => Some(Element::Er),
            "Tm" => Some(Element::Tm),
            "Yb" => Some(Element::Yb),
            "Lu" => Some(Element::Lu),
            "Hf" => Some(Element::Hf),
            "Ta" => Some(Element::Ta),
            "W" => Some(Element::W),
            "Re" => Some(Element::Re),
            "Os" => Some(Element::Os),
            "Ir" => Some(Element::Ir),
            "Pt" => Some(Element::Pt),
            "Au" => Some(Element::Au),
            "Hg" => Some(Element::Hg),
            "Tl" => Some(Element::Tl),
            "Pb" => Some(Element::Pb),
            "Bi" => Some(Element::Bi),
            "Po" => Some(Element::Po),
            "At" => Some(Element::At),
            "Rn" => Some(Element::Rn),
            "Fr" => Some(Element::Fr),
            "Ra" => Some(Element::Ra),
            "Ac" => Some(Element::Ac),
            "Th" => Some(Element::Th),
            "Pa" => Some(Element::Pa),
            "U" => Some(Element::U),
            "Np" => Some(Element::Np),
            "Pu" => Some(Element::Pu),
            "Am" => Some(Element::Am),
            "Cm" => Some(Element::Cm),
            "Bk" => Some(Element::Bk),
            "Cf" => Some(Element::Cf),
            "Es" => Some(Element::Es),
            "Fm" => Some(Element::Fm),
            "Md" => Some(Element::Md),
            "No" => Some(Element::No),
            "Lr" => Some(Element::Lr),
            "Rf" => Some(Element::Rf),
            "Db" => Some(Element::Db),
            "Sg" => Some(Element::Sg),
            "Bh" => Some(Element::Bh),
            "Hs" => Some(Element::Hs),
            "Mt" => Some(Element::Mt),
            "Ds" => Some(Element::Ds),
            "Rg" => Some(Element::Rg),
            "Cn" => Some(Element::Cn),
            "Nh" => Some(Element::Nh),
            "Fl" => Some(Element::Fl),
            "Mc" => Some(Element::Mc),
            "Lv" => Some(Element::Lv),
            "Ts" => Some(Element::Ts),
            "Og" => Some(Element::Og),
            _ => None,
        }
    }

    pub const fn properties(self) -> Properties {
        match self {
            Element::H => {
                Properties { symbol: "H", name: "Hydrogen", weight: 1.008, covalent_radius: 0.31, van_der_waals_radius: Some(1.20), valences: &[1] }
            }
            Element::He => {
                Properties { symbol: "He", name: "Helium", weight: 4.0026, covalent_radius: 0.28, van_der_waals_radius: Some(1.43), valences: &[0] }
            }
            Element::Li => {
                Properties { symbol: "Li", name: "Lithium", weight: 6.94, covalent_radius: 1.28, van_der_waals_radius: Some(2.12), valences: &[] }
            }
            Element::Be => Properties {
                symbol: "Be",
                name: "Beryllium",
                weight: 9.0122,
                covalent_radius: 0.96,
                van_der_waals_radius: Some(1.98),
                valences: &[2],
            },
            Element::B => {
                Properties { symbol: "B", name: "Boron", weight: 10.81, covalent_radius: 0.84, van_der_waals_radius: Some(1.91), valences: &[3] }
            }
            Element::C => {
                Properties { symbol: "C", name: "Carbon", weight: 12.011, covalent_radius: 0.76, van_der_waals_radius: Some(1.77), valences: &[4] }
            }
            Element::N => {
                Properties { symbol: "N", name: "Nitrogen", weight: 14.007, covalent_radius: 0.71, van_der_waals_radius: Some(1.66), valences: &[3] }
            }
            Element::O => {
                Properties { symbol: "O", name: "Oxygen", weight: 15.999, covalent_radius: 0.66, van_der_waals_radius: Some(1.50), valences: &[2] }
            }
            Element::F => {
                Properties { symbol: "F", name: "Fluorine", weight: 18.998, covalent_radius: 0.57, van_der_waals_radius: Some(1.46), valences: &[1] }
            }
            Element::Ne => {
                Properties { symbol: "Ne", name: "Neon", weight: 20.180, covalent_radius: 0.58, van_der_waals_radius: Some(1.58), valences: &[0] }
            }
            Element::Na => {
                Properties { symbol: "Na", name: "Sodium", weight: 22.990, covalent_radius: 1.66, van_der_waals_radius: Some(2.50), valences: &[] }
            }
            Element::Mg => {
                Properties { symbol: "Mg", name: "Magnesium", weight: 24.305, covalent_radius: 1.41, van_der_waals_radius: Some(2.51), valences: &[] }
            }
            Element::Al => Properties {
                symbol: "Al",
                name: "Aluminium",
                weight: 26.982,
                covalent_radius: 1.21,
                van_der_waals_radius: Some(2.25),
                valences: &[3],
            },
            Element::Si => {
                Properties { symbol: "Si", name: "Silicon", weight: 28.085, covalent_radius: 1.11, van_der_waals_radius: Some(2.19), valences: &[4] }
            }
            Element::P => Properties {
                symbol: "P",
                name: "Phosphorus",
                weight: 30.974,
                covalent_radius: 1.07,
                van_der_waals_radius: Some(1.90),
                valences: &[3, 5],
            },
            Element::S => Properties {
                symbol: "S",
                name: "Sulfur",
                weight: 32.06,
                covalent_radius: 1.05,
                van_der_waals_radius: Some(1.89),
                valences: &[2, 4, 6],
            },
            Element::Cl => {
                Properties { symbol: "Cl", name: "Chlorine", weight: 35.45, covalent_radius: 1.02, van_der_waals_radius: Some(1.82), valences: &[1] }
            }
            Element::Ar => {
                Properties { symbol: "Ar", name: "Argon", weight: 39.95, covalent_radius: 1.06, van_der_waals_radius: Some(1.94), valences: &[0] }
            }
            Element::K => {
                Properties { symbol: "K", name: "Potassium", weight: 39.098, covalent_radius: 2.03, van_der_waals_radius: Some(2.73), valences: &[] }
            }
            Element::Ca => {
                Properties { symbol: "Ca", name: "Calcium", weight: 40.078, covalent_radius: 1.76, van_der_waals_radius: Some(2.62), valences: &[] }
            }
            Element::Sc => {
                Properties { symbol: "Sc", name: "Scandium", weight: 44.956, covalent_radius: 1.70, van_der_waals_radius: Some(2.58), valences: &[] }
            }
            Element::Ti => {
                Properties { symbol: "Ti", name: "Titanium", weight: 47.867, covalent_radius: 1.60, van_der_waals_radius: Some(2.46), valences: &[] }
            }
            Element::V => {
                Properties { symbol: "V", name: "Vanadium", weight: 50.942, covalent_radius: 1.53, van_der_waals_radius: Some(2.42), valences: &[] }
            }
            Element::Cr => {
                Properties { symbol: "Cr", name: "Chromium", weight: 51.996, covalent_radius: 1.39, van_der_waals_radius: Some(2.45), valences: &[] }
            }
            Element::Mn => {
                Properties { symbol: "Mn", name: "Manganese", weight: 54.938, covalent_radius: 1.39, van_der_waals_radius: Some(2.45), valences: &[] }
            }
            Element::Fe => {
                Properties { symbol: "Fe", name: "Iron", weight: 55.845, covalent_radius: 1.32, van_der_waals_radius: Some(2.44), valences: &[] }
            }
            Element::Co => {
                Properties { symbol: "Co", name: "Cobalt", weight: 58.933, covalent_radius: 1.26, van_der_waals_radius: Some(2.40), valences: &[] }
            }
            Element::Ni => {
                Properties { symbol: "Ni", name: "Nickel", weight: 58.693, covalent_radius: 1.24, van_der_waals_radius: Some(2.40), valences: &[] }
            }
            Element::Cu => {
                Properties { symbol: "Cu", name: "Copper", weight: 63.546, covalent_radius: 1.32, van_der_waals_radius: Some(2.38), valences: &[] }
            }
            Element::Zn => {
                Properties { symbol: "Zn", name: "Zinc", weight: 65.38, covalent_radius: 1.22, van_der_waals_radius: Some(2.39), valences: &[] }
            }
            Element::Ga => {
                Properties { symbol: "Ga", name: "Gallium", weight: 69.723, covalent_radius: 1.22, van_der_waals_radius: Some(2.32), valences: &[3] }
            }
            Element::Ge => Properties {
                symbol: "Ge",
                name: "Germanium",
                weight: 72.630,
                covalent_radius: 1.20,
                van_der_waals_radius: Some(2.29),
                valences: &[4],
            },
            Element::As => Properties {
                symbol: "As",
                name: "Arsenic",
                weight: 74.922,
                covalent_radius: 1.19,
                van_der_waals_radius: Some(1.88),
                valences: &[3, 5],
            },
            Element::Se => Properties {
                symbol: "Se",
                name: "Selenium",
                weight: 78.971,
                covalent_radius: 1.20,
                van_der_waals_radius: Some(1.82),
                valences: &[2, 4, 6],
            },
            Element::Br => {
                Properties { symbol: "Br", name: "Bromine", weight: 79.904, covalent_radius: 1.20, van_der_waals_radius: Some(1.86), valences: &[1] }
            }
            Element::Kr => {
                Properties { symbol: "Kr", name: "Krypton", weight: 83.798, covalent_radius: 1.16, van_der_waals_radius: Some(2.07), valences: &[0] }
            }
            Element::Rb => {
                Properties { symbol: "Rb", name: "Rubidium", weight: 85.468, covalent_radius: 2.20, van_der_waals_radius: Some(3.21), valences: &[] }
            }
            Element::Sr => {
                Properties { symbol: "Sr", name: "Strontium", weight: 87.62, covalent_radius: 1.95, van_der_waals_radius: Some(2.84), valences: &[] }
            }
            Element::Y => {
                Properties { symbol: "Y", name: "Yttrium", weight: 88.906, covalent_radius: 1.90, van_der_waals_radius: Some(2.75), valences: &[] }
            }
            Element::Zr => {
                Properties { symbol: "Zr", name: "Zirconium", weight: 91.224, covalent_radius: 1.75, van_der_waals_radius: Some(2.52), valences: &[] }
            }
            Element::Nb => {
                Properties { symbol: "Nb", name: "Niobium", weight: 92.906, covalent_radius: 1.64, van_der_waals_radius: Some(2.56), valences: &[] }
            }
            Element::Mo => {
                Properties { symbol: "Mo", name: "Molybdenum", weight: 95.95, covalent_radius: 1.54, van_der_waals_radius: Some(2.45), valences: &[] }
            }
            Element::Tc => {
                Properties { symbol: "Tc", name: "Technetium", weight: 98.0, covalent_radius: 1.47, van_der_waals_radius: Some(2.44), valences: &[] }
            }
            Element::Ru => {
                Properties { symbol: "Ru", name: "Ruthenium", weight: 101.07, covalent_radius: 1.46, van_der_waals_radius: Some(2.46), valences: &[] }
            }
            Element::Rh => {
                Properties { symbol: "Rh", name: "Rhodium", weight: 102.91, covalent_radius: 1.42, van_der_waals_radius: Some(2.44), valences: &[] }
            }
            Element::Pd => {
                Properties { symbol: "Pd", name: "Palladium", weight: 106.42, covalent_radius: 1.39, van_der_waals_radius: Some(2.15), valences: &[] }
            }
            Element::Ag => {
                Properties { symbol: "Ag", name: "Silver", weight: 107.87, covalent_radius: 1.45, van_der_waals_radius: Some(2.53), valences: &[] }
            }
            Element::Cd => {
                Properties { symbol: "Cd", name: "Cadmium", weight: 112.41, covalent_radius: 1.44, van_der_waals_radius: Some(2.49), valences: &[] }
            }
            Element::In => {
                Properties { symbol: "In", name: "Indium", weight: 114.82, covalent_radius: 1.42, van_der_waals_radius: Some(2.43), valences: &[3] }
            }
            Element::Sn => {
                Properties { symbol: "Sn", name: "Tin", weight: 118.71, covalent_radius: 1.39, van_der_waals_radius: Some(2.42), valences: &[2, 4] }
            }
            Element::Sb => Properties {
                symbol: "Sb",
                name: "Antimony",
                weight: 121.76,
                covalent_radius: 1.39,
                van_der_waals_radius: Some(2.47),
                valences: &[3, 5],
            },
            Element::Te => Properties {
                symbol: "Te",
                name: "Tellurium",
                weight: 127.60,
                covalent_radius: 1.38,
                van_der_waals_radius: Some(1.99),
                valences: &[2, 4, 6],
            },
            Element::I => Properties {
                symbol: "I",
                name: "Iodine",
                weight: 126.90,
                covalent_radius: 1.39,
                van_der_waals_radius: Some(2.04),
                valences: &[1, 3, 5],
            },
            Element::Xe => Properties {
                symbol: "Xe",
                name: "Xenon",
                weight: 131.29,
                covalent_radius: 1.40,
                van_der_waals_radius: Some(2.28),
                valences: &[0, 2, 4, 6],
            },
            Element::Cs => {
                Properties { symbol: "Cs", name: "Caesium", weight: 132.91, covalent_radius: 2.44, van_der_waals_radius: Some(3.48), valences: &[1] }
            }
            Element::Ba => {
                Properties { symbol: "Ba", name: "Barium", weight: 137.33, covalent_radius: 2.15, van_der_waals_radius: Some(3.03), valences: &[] }
            }
            Element::La => {
                Properties { symbol: "La", name: "Lanthanum", weight: 138.91, covalent_radius: 2.07, van_der_waals_radius: Some(2.98), valences: &[] }
            }
            Element::Ce => {
                Properties { symbol: "Ce", name: "Cerium", weight: 140.12, covalent_radius: 2.04, van_der_waals_radius: Some(2.88), valences: &[] }
            }
            Element::Pr => Properties {
                symbol: "Pr",
                name: "Praseodymium",
                weight: 140.91,
                covalent_radius: 2.03,
                van_der_waals_radius: Some(2.92),
                valences: &[],
            },
            Element::Nd => {
                Properties { symbol: "Nd", name: "Neodymium", weight: 144.24, covalent_radius: 2.01, van_der_waals_radius: Some(2.95), valences: &[] }
            }
            Element::Pm => {
                Properties { symbol: "Pm", name: "Promethium", weight: 145.0, covalent_radius: 1.99, van_der_waals_radius: None, valences: &[] }
            }
            Element::Sm => {
                Properties { symbol: "Sm", name: "Samarium", weight: 150.36, covalent_radius: 1.98, van_der_waals_radius: Some(2.90), valences: &[] }
            }
            Element::Eu => {
                Properties { symbol: "Eu", name: "Europium", weight: 151.96, covalent_radius: 1.98, van_der_waals_radius: Some(2.87), valences: &[] }
            }
            Element::Gd => Properties {
                symbol: "Gd",
                name: "Gadolinium",
                weight: 157.25,
                covalent_radius: 1.96,
                van_der_waals_radius: Some(2.83),
                valences: &[],
            },
            Element::Tb => {
                Properties { symbol: "Tb", name: "Terbium", weight: 158.93, covalent_radius: 1.94, van_der_waals_radius: Some(2.79), valences: &[] }
            }
            Element::Dy => Properties {
                symbol: "Dy",
                name: "Dysprosium",
                weight: 162.50,
                covalent_radius: 1.92,
                van_der_waals_radius: Some(2.87),
                valences: &[],
            },
            Element::Ho => {
                Properties { symbol: "Ho", name: "Holmium", weight: 164.93, covalent_radius: 1.92, van_der_waals_radius: Some(2.81), valences: &[] }
            }
            Element::Er => {
                Properties { symbol: "Er", name: "Erbium", weight: 167.26, covalent_radius: 1.89, van_der_waals_radius: Some(2.83), valences: &[] }
            }
            Element::Tm => {
                Properties { symbol: "Tm", name: "Thulium", weight: 168.93, covalent_radius: 1.90, van_der_waals_radius: Some(2.79), valences: &[] }
            }
            Element::Yb => {
                Properties { symbol: "Yb", name: "Ytterbium", weight: 173.05, covalent_radius: 1.87, van_der_waals_radius: Some(2.80), valences: &[] }
            }
            Element::Lu => {
                Properties { symbol: "Lu", name: "Lutetium", weight: 174.97, covalent_radius: 1.87, van_der_waals_radius: Some(2.74), valences: &[] }
            }
            Element::Hf => {
                Properties { symbol: "Hf", name: "Hafnium", weight: 178.49, covalent_radius: 1.75, van_der_waals_radius: Some(2.63), valences: &[] }
            }
            Element::Ta => {
                Properties { symbol: "Ta", name: "Tantalum", weight: 180.95, covalent_radius: 1.70, van_der_waals_radius: Some(2.53), valences: &[] }
            }
            Element::W => {
                Properties { symbol: "W", name: "Tungsten", weight: 183.84, covalent_radius: 1.62, van_der_waals_radius: Some(2.57), valences: &[] }
            }
            Element::Re => {
                Properties { symbol: "Re", name: "Rhenium", weight: 186.21, covalent_radius: 1.51, van_der_waals_radius: Some(2.49), valences: &[] }
            }
            Element::Os => {
                Properties { symbol: "Os", name: "Osmium", weight: 190.23, covalent_radius: 1.44, van_der_waals_radius: Some(2.48), valences: &[] }
            }
            Element::Ir => {
                Properties { symbol: "Ir", name: "Iridium", weight: 192.22, covalent_radius: 1.41, van_der_waals_radius: Some(2.41), valences: &[] }
            }
            Element::Pt => {
                Properties { symbol: "Pt", name: "Platinum", weight: 195.08, covalent_radius: 1.36, van_der_waals_radius: Some(2.29), valences: &[] }
            }
            Element::Au => {
                Properties { symbol: "Au", name: "Gold", weight: 196.97, covalent_radius: 1.36, van_der_waals_radius: Some(2.32), valences: &[] }
            }
            Element::Hg => {
                Properties { symbol: "Hg", name: "Mercury", weight: 200.59, covalent_radius: 1.32, van_der_waals_radius: Some(2.45), valences: &[] }
            }
            Element::Tl => {
                Properties { symbol: "Tl", name: "Thallium", weight: 204.38, covalent_radius: 1.45, van_der_waals_radius: Some(2.47), valences: &[] }
            }
            Element::Pb => {
                Properties { symbol: "Pb", name: "Lead", weight: 207.2, covalent_radius: 1.46, van_der_waals_radius: Some(2.60), valences: &[2, 4] }
            }
            Element::Bi => Properties {
                symbol: "Bi",
                name: "Bismuth",
                weight: 208.98,
                covalent_radius: 1.48,
                van_der_waals_radius: Some(2.54),
                valences: &[3, 5],
            },
            Element::Po => {
                Properties { symbol: "Po", name: "Polonium", weight: 209.0, covalent_radius: 1.40, van_der_waals_radius: None, valences: &[2, 4, 6] }
            }
            Element::At => {
                Properties { symbol: "At", name: "Astatine", weight: 210.0, covalent_radius: 1.50, van_der_waals_radius: None, valences: &[1, 3, 5] }
            }
            Element::Rn => {
                Properties { symbol: "Rn", name: "Radon", weight: 222.0, covalent_radius: 1.50, van_der_waals_radius: Some(2.40), valences: &[0] }
            }
            Element::Fr => {
                Properties { symbol: "Fr", name: "Francium", weight: 223.0, covalent_radius: 2.60, van_der_waals_radius: None, valences: &[1] }
            }
            Element::Ra => {
                Properties { symbol: "Ra", name: "Radium", weight: 226.0, covalent_radius: 2.21, van_der_waals_radius: None, valences: &[] }
            }
            Element::Ac => {
                Properties { symbol: "Ac", name: "Actinium", weight: 227.0, covalent_radius: 2.15, van_der_waals_radius: Some(2.80), valences: &[] }
            }
            Element::Th => {
                Properties { symbol: "Th", name: "Thorium", weight: 232.04, covalent_radius: 2.06, van_der_waals_radius: Some(2.93), valences: &[] }
            }
            Element::Pa => Properties {
                symbol: "Pa",
                name: "Protactinium",
                weight: 231.04,
                covalent_radius: 2.00,
                van_der_waals_radius: Some(2.88),
                valences: &[],
            },
            Element::U => {
                Properties { symbol: "U", name: "Uranium", weight: 238.03, covalent_radius: 1.96, van_der_waals_radius: Some(2.71), valences: &[] }
            }
            Element::Np => {
                Properties { symbol: "Np", name: "Neptunium", weight: 237.0, covalent_radius: 1.90, van_der_waals_radius: Some(2.82), valences: &[] }
            }
            Element::Pu => {
                Properties { symbol: "Pu", name: "Plutonium", weight: 244.0, covalent_radius: 1.87, van_der_waals_radius: Some(2.81), valences: &[] }
            }
            Element::Am => {
                Properties { symbol: "Am", name: "Americium", weight: 243.0, covalent_radius: 1.80, van_der_waals_radius: Some(2.83), valences: &[] }
            }
            Element::Cm => {
                Properties { symbol: "Cm", name: "Curium", weight: 247.0, covalent_radius: 1.69, van_der_waals_radius: Some(3.05), valences: &[] }
            }
            Element::Bk => {
                Properties { symbol: "Bk", name: "Berkelium", weight: 247.0, covalent_radius: 1.68, van_der_waals_radius: Some(3.40), valences: &[] }
            }
            Element::Cf => Properties {
                symbol: "Cf",
                name: "Californium",
                weight: 251.0,
                covalent_radius: 1.68,
                van_der_waals_radius: Some(3.05),
                valences: &[],
            },
            Element::Es => Properties {
                symbol: "Es",
                name: "Einsteinium",
                weight: 252.0,
                covalent_radius: 1.65,
                van_der_waals_radius: Some(2.70),
                valences: &[],
            },
            Element::Fm => {
                Properties { symbol: "Fm", name: "Fermium", weight: 257.0, covalent_radius: 1.67, van_der_waals_radius: None, valences: &[] }
            }
            Element::Md => {
                Properties { symbol: "Md", name: "Mendelevium", weight: 258.0, covalent_radius: 1.73, van_der_waals_radius: None, valences: &[] }
            }
            Element::No => {
                Properties { symbol: "No", name: "Nobelium", weight: 259.0, covalent_radius: 1.76, van_der_waals_radius: None, valences: &[] }
            }
            Element::Lr => {
                Properties { symbol: "Lr", name: "Lawrencium", weight: 266.0, covalent_radius: 1.61, van_der_waals_radius: None, valences: &[] }
            }
            Element::Rf => {
                Properties { symbol: "Rf", name: "Rutherfordium", weight: 267.0, covalent_radius: 1.57, van_der_waals_radius: None, valences: &[] }
            }
            Element::Db => {
                Properties { symbol: "Db", name: "Dubnium", weight: 268.0, covalent_radius: 1.49, van_der_waals_radius: None, valences: &[] }
            }
            Element::Sg => {
                Properties { symbol: "Sg", name: "Seaborgium", weight: 269.0, covalent_radius: 1.43, van_der_waals_radius: None, valences: &[] }
            }
            Element::Bh => {
                Properties { symbol: "Bh", name: "Bohrium", weight: 270.0, covalent_radius: 1.41, van_der_waals_radius: None, valences: &[] }
            }
            Element::Hs => {
                Properties { symbol: "Hs", name: "Hassium", weight: 269.0, covalent_radius: 1.34, van_der_waals_radius: None, valences: &[] }
            }
            Element::Mt => {
                Properties { symbol: "Mt", name: "Meitnerium", weight: 278.0, covalent_radius: 1.29, van_der_waals_radius: None, valences: &[] }
            }
            Element::Ds => {
                Properties { symbol: "Ds", name: "Darmstadtium", weight: 281.0, covalent_radius: 1.28, van_der_waals_radius: None, valences: &[] }
            }
            Element::Rg => {
                Properties { symbol: "Rg", name: "Roentgenium", weight: 282.0, covalent_radius: 1.21, van_der_waals_radius: None, valences: &[] }
            }
            Element::Cn => {
                Properties { symbol: "Cn", name: "Copernicium", weight: 285.0, covalent_radius: 1.22, van_der_waals_radius: None, valences: &[] }
            }
            Element::Nh => {
                Properties { symbol: "Nh", name: "Nihonium", weight: 286.0, covalent_radius: 1.36, van_der_waals_radius: None, valences: &[] }
            }
            Element::Fl => {
                Properties { symbol: "Fl", name: "Flerovium", weight: 289.0, covalent_radius: 1.43, van_der_waals_radius: None, valences: &[] }
            }
            Element::Mc => {
                Properties { symbol: "Mc", name: "Moscovium", weight: 290.0, covalent_radius: 1.62, van_der_waals_radius: None, valences: &[] }
            }
            Element::Lv => {
                Properties { symbol: "Lv", name: "Livermorium", weight: 293.0, covalent_radius: 1.75, van_der_waals_radius: None, valences: &[] }
            }
            Element::Ts => {
                Properties { symbol: "Ts", name: "Tennessine", weight: 294.0, covalent_radius: 1.65, van_der_waals_radius: None, valences: &[] }
            }
            Element::Og => {
                Properties { symbol: "Og", name: "Oganesson", weight: 294.0, covalent_radius: 1.57, van_der_waals_radius: None, valences: &[] }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn discriminant_equals_atomic_number() {
        assert_eq!(Element::H.atomic_number(), 1);
        assert_eq!(Element::C.atomic_number(), 6);
        assert_eq!(Element::Ca.atomic_number(), 20);
        assert_eq!(Element::Og.atomic_number(), 118);
    }

    #[test]
    fn properties_match_reports_own_symbol() {
        // The `properties()` match is keyed by variant; this catches a
        // copy-paste where an arm returns the wrong element's data.
        assert_eq!(Element::N.symbol(), "N");
        assert_eq!(Element::N.weight(), 14.007);
        assert_eq!(Element::Ca.symbol(), "Ca");
    }

    #[test]
    fn from_symbol_round_trips_with_symbol() {
        assert_eq!(Element::from_symbol("C"), Some(Element::C));
        assert_eq!(Element::from_symbol("Ca"), Some(Element::Ca));
        assert_eq!(Element::from_symbol("Og"), Some(Element::Og));
        assert_eq!(Element::from_symbol("Xx"), None);
    }

    #[test]
    fn van_der_waals_radius_is_larger_than_covalent_radius() {
        assert_eq!(Element::H.van_der_waals_radius(), Some(1.20));
        assert_eq!(Element::Au.van_der_waals_radius(), Some(2.32));
        assert!(Element::C.van_der_waals_radius().unwrap() > Element::C.covalent_radius());
    }

    #[test]
    fn valences_list_every_normal_bonding_state() {
        assert_eq!(Element::C.valences(), &[4]);
        assert_eq!(Element::N.valences(), &[3]);
        assert_eq!(Element::P.valences(), &[3, 5]);
        assert_eq!(Element::S.valences(), &[2, 4, 6]);
        assert_eq!(Element::Ne.valences(), &[0]);
    }

    #[test]
    fn valences_are_unconstrained_for_metals() {
        assert!(Element::Ca.valences().is_empty());
        assert!(Element::Fe.valences().is_empty());
        assert!(Element::Zn.valences().is_empty());
    }

    #[test]
    fn van_der_waals_radius_is_none_for_undefined_elements() {
        assert_eq!(Element::Pm.van_der_waals_radius(), None);
        assert_eq!(Element::Og.van_der_waals_radius(), None);
    }
}
