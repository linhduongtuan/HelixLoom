//! Conversions between mass, moles, concentration and volume, with the molecular weight of the substance (typed, or worked
//! out from a DNA or RNA sequence). Re-written from `calc_molarity` (AppMain.tcl:35237, 0x19aec7) and `calc_set_units`
//! (L35463, 0x19dae4) with the arithmetic in the same order.

use super::format::{g3, number};

/// A unit of the first and second quantity (`g mg ug ng pg mole mmole umole nmole pmole`, and `%(wt/vol)` for the first).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Unit {
    G,
    Mg,
    Ug,
    Ng,
    Pg,
    Mole,
    Mmole,
    Umole,
    Nmole,
    Pmole,
    PercentWv,
}

impl Unit {
    /// The menu of the first quantity.
    pub const FIRST: [Unit; 11] = [
        Unit::G,
        Unit::Mg,
        Unit::Ug,
        Unit::Ng,
        Unit::Pg,
        Unit::Mole,
        Unit::Mmole,
        Unit::Umole,
        Unit::Nmole,
        Unit::Pmole,
        Unit::PercentWv,
    ];
    /// The menu of the second quantity: the same without the percentage.
    pub const SECOND: [Unit; 10] = [
        Unit::G,
        Unit::Mg,
        Unit::Ug,
        Unit::Ng,
        Unit::Pg,
        Unit::Mole,
        Unit::Mmole,
        Unit::Umole,
        Unit::Nmole,
        Unit::Pmole,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Unit::G => "g",
            Unit::Mg => "mg",
            Unit::Ug => "ug",
            Unit::Ng => "ng",
            Unit::Pg => "pg",
            Unit::Mole => "mole",
            Unit::Mmole => "mmole",
            Unit::Umole => "umole",
            Unit::Nmole => "nmole",
            Unit::Pmole => "pmole",
            Unit::PercentWv => "%(wt/vol)",
        }
    }

    /// How many grams or moles one of it is (`unit_conv1`, `unit_conv2`).
    fn factor(self) -> f64 {
        match self {
            Unit::PercentWv | Unit::G | Unit::Mole => 1.0,
            Unit::Mg | Unit::Mmole => 1e-3,
            Unit::Ug | Unit::Umole => 1e-6,
            Unit::Ng | Unit::Nmole => 1e-9,
            Unit::Pg | Unit::Pmole => 1e-12,
        }
    }

    fn is_mass(self) -> bool {
        matches!(self, Unit::G | Unit::Mg | Unit::Ug | Unit::Ng | Unit::Pg)
    }

    fn is_mole(self) -> bool {
        matches!(
            self,
            Unit::Mole | Unit::Mmole | Unit::Umole | Unit::Nmole | Unit::Pmole
        )
    }
}

/// The volume a quantity is per (`/L /mL /uL`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Per {
    L,
    Ml,
    Ul,
}

impl Per {
    pub const ALL: [Per; 3] = [Per::L, Per::Ml, Per::Ul];

    pub fn label(self) -> &'static str {
        match self {
            Per::L => "/L",
            Per::Ml => "/mL",
            Per::Ul => "/uL",
        }
    }

    fn factor(self) -> f64 {
        match self {
            Per::L => 1.0,
            Per::Ml => 1e3,
            Per::Ul => 1e6,
        }
    }
}

/// Where the molecular weight comes from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Weight {
    /// Typed, in g/mol.
    GramsPerMole,
    /// Single-stranded DNA: typed as a length, or read from a sequence.
    SsDna,
    DsDna,
    Rna,
}

impl Weight {
    pub const ALL: [Weight; 4] = [
        Weight::GramsPerMole,
        Weight::SsDna,
        Weight::DsDna,
        Weight::Rna,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Weight::GramsPerMole => "g/mol",
            Weight::SsDna => "nucleotides ssDNA",
            Weight::DsDna => "base pairs dsDNA",
            Weight::Rna => "nucleotides RNA",
        }
    }

    /// A sequence can be typed for these.
    pub fn takes_sequence(self) -> bool {
        self != Weight::GramsPerMole
    }
}

/// What the second line gives its amount in (`L mL uL mole mmole umole nmole g mg ug ng`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Amount {
    L,
    Ml,
    Ul,
    Mole,
    Mmole,
    Umole,
    Nmole,
    G,
    Mg,
    Ug,
    Ng,
}

/// Whether an amount is a volume, moles or grams.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Volume,
    Mole,
    Mass,
}

impl Amount {
    pub const ALL: [Amount; 11] = [
        Amount::L,
        Amount::Ml,
        Amount::Ul,
        Amount::Mole,
        Amount::Mmole,
        Amount::Umole,
        Amount::Nmole,
        Amount::G,
        Amount::Mg,
        Amount::Ug,
        Amount::Ng,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Amount::L => "L",
            Amount::Ml => "mL",
            Amount::Ul => "uL",
            Amount::Mole => "mole",
            Amount::Mmole => "mmole",
            Amount::Umole => "umole",
            Amount::Nmole => "nmole",
            Amount::G => "g",
            Amount::Mg => "mg",
            Amount::Ug => "ug",
            Amount::Ng => "ng",
        }
    }

    /// How many of it make one litre, mole or gram (`unit_conv3`).
    fn factor(self) -> f64 {
        match self {
            Amount::L | Amount::Mole | Amount::G => 1.0,
            Amount::Ml | Amount::Mmole | Amount::Mg => 1e3,
            Amount::Ul | Amount::Umole | Amount::Ug => 1e6,
            Amount::Nmole | Amount::Ng => 1e9,
        }
    }

    fn kind(self) -> Kind {
        match self {
            Amount::L | Amount::Ml | Amount::Ul => Kind::Volume,
            Amount::Mole | Amount::Mmole | Amount::Umole | Amount::Nmole => Kind::Mole,
            Amount::G | Amount::Mg | Amount::Ug | Amount::Ng => Kind::Mass,
        }
    }
}

/// The weights of the nucleotides in g/mol (each with a water off, the ends add the water back): DNA singles A 313.23,
/// C 289.20, G 329.23, T 304.21 (+18.02); base pairs AT 617.44, GC 618.22 (+36.04); RNA A 329.2, C 305.2, G 345.2,
/// U 306.2 (+18.02). Without a sequence an average is used: 308.97, 617.96 and 321.47 for each one.
const SS_DNA: (f64, f64, f64, f64, f64) = (313.23, 289.20, 329.23, 304.21, 18.02);
const DS_DNA: (f64, f64, f64) = (617.44, 618.22, 36.04);
const RNA: (f64, f64, f64, f64, f64) = (329.2, 305.2, 345.2, 306.2, 18.02);
const SS_DNA_AVERAGE: f64 = 308.97;
const DS_DNA_AVERAGE: f64 = 617.96;
const RNA_AVERAGE: f64 = 321.47;

#[derive(Clone, Copy, Debug)]
pub struct MolarityInput<'a> {
    pub first: &'a str,
    pub unit1: Unit,
    /// The volume the first quantity is per: none for a plain amount. Ignored for a percentage.
    pub per1: Option<Per>,
    pub unit2: Unit,
    pub per2: Per,
    pub weight: Weight,
    /// The molecular weight in g/mol, or the length for the DNA and RNA ones.
    pub mol_wt: &'a str,
    pub sequence: &'a str,
    /// The amount the second line is worked out for, and what it is given in.
    pub amount: &'a str,
    pub amount_unit: Amount,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MolarityOutput {
    /// The first quantity in the second unit (`%.3g`); empty when it cannot be worked out.
    pub second: String,
    /// The molecular weight field after a sequence was read: its length.
    pub mol_wt: String,
    /// The sequence as it was kept: only A, C, G and T (U for RNA).
    pub sequence: String,
    /// The second line: what the amount comes to; empty when there is none to show.
    pub detail: String,
    /// Which kind of second line it is: 3, 5 or 7 (moles, grams, or both), 0 for none.
    pub key: u8,
    /// A molecular weight is needed to convert between the two units.
    pub needs_weight: bool,
    /// The quantity is per a volume, so the second line is of use.
    pub per_volume: bool,
}

/// `format %.3g` of a value with the metric prefix that suits it (`calc_set_units`): `1 mg`, `0.1 pmole`.
pub fn set_units(value: f64, units: &str) -> Option<String> {
    if value.is_nan() || value < 0.0 {
        return None;
    }
    // the decimal exponent is found with natural logs as the original does, so 1e-3 is 3 below zero and a bit
    let l = (value.ln() / 10f64.ln()).floor();
    Some(if l < -9.0 {
        format!("{} p{units}", g3(value * 1e12))
    } else if l < -6.0 {
        format!("{} n{units}", g3(value * 1e9))
    } else if l < -3.0 {
        format!("{} u{units}", g3(value * 1e6))
    } else if l < 0.0 {
        format!("{} m{units}", g3(value * 1e3))
    } else {
        format!("{} {units}", g3(value))
    })
}

/// The letters kept of a sequence: U is T for DNA and T is U for RNA, and nothing but A, C, G and T (U) is left.
fn clean_sequence(text: &str, rna: bool) -> String {
    text.chars()
        .map(|c| match (rna, c) {
            (true, 'T') => 'U',
            (true, 't') => 'u',
            (false, 'U') => 'T',
            (false, 'u') => 't',
            _ => c,
        })
        .filter(|c| "ACGTUacgtu".contains(*c))
        .collect()
}

/// How many of each letter (A, C, G, T or U) a cleaned sequence holds, in any case.
fn counts(sequence: &str) -> (f64, f64, f64, f64) {
    let count = |set: &str| sequence.chars().filter(|c| set.contains(*c)).count() as f64;
    (count("Aa"), count("Cc"), count("Gg"), count("TtUu"))
}

fn molecular_weight(weight: Weight, sequence: &str, length: f64) -> f64 {
    let (a, c, g, t) = counts(sequence);
    let has_sequence = !sequence.is_empty();
    match weight {
        Weight::GramsPerMole => 1.0 * length,
        Weight::SsDna if has_sequence => {
            let (wa, wc, wg, wt, water) = SS_DNA;
            1.0 * a * wa + c * wc + g * wg + t * wt + water
        }
        Weight::SsDna => 1.0 * length * SS_DNA_AVERAGE + SS_DNA.4,
        Weight::DsDna if has_sequence => {
            let (at_pair, gc_pair, ends) = DS_DNA;
            1.0 * (a + t) * at_pair + (c + g) * gc_pair + ends
        }
        Weight::DsDna => 1.0 * length * DS_DNA_AVERAGE + DS_DNA.2,
        Weight::Rna if has_sequence => {
            let (wa, wc, wg, wu, water) = RNA;
            1.0 * a * wa + c * wc + g * wg + t * wu + water
        }
        Weight::Rna => 1.0 * length * RNA_AVERAGE + RNA.4,
    }
}

/// The second line of the original: what the amount comes to, in grams and moles or litres, given the first quantity
/// worked out as a number of grams or moles per litre (`calc_value`).
///
/// **Deliberate difference:** with moles per litre and no mass (the original's case 3) the original divides by the
/// volume where every other case multiplies (1 mmole/L in 0.5 L is "2 mmole"), and its other two branches test the
/// amount instead of its unit, so they never run. Here the amount is the concentration times the volume, and the volume
/// of an amount of moles is the moles over the concentration.
fn second_line(
    key: u8,
    amount_unit: Amount,
    volume: f64,
    calc_value: f64,
    weight: f64,
) -> Option<String> {
    let litres = volume / amount_unit.factor();
    let kind = amount_unit.kind();
    match (key, kind) {
        (3, Kind::Volume) => Some(format!(
            " =  {} ",
            set_units(1.0 * calc_value * litres, "mole")?
        )),
        (3, Kind::Mole) => Some(format!(
            " =  {} ",
            set_units(1.0 * litres / calc_value, "L")?
        )),
        (5, Kind::Volume) => Some(format!(
            " =  {} ",
            set_units(1.0 * calc_value * litres, "g")?
        )),
        (5, Kind::Mass) => Some(format!(" = {}", set_units(1.0 * litres / calc_value, "L")?)),
        (5, Kind::Mole) => Some("   ".into()),
        (7, Kind::Volume) => Some(format!(
            " =  {} , {} ",
            set_units(1.0 * calc_value * litres, "g")?,
            set_units(1.0 * calc_value * litres / weight, "mole")?
        )),
        (7, Kind::Mass) => Some(format!(
            " = {}, {}",
            set_units(1.0 * litres / calc_value, "L")?,
            set_units(1.0 * litres / weight, "mole")?
        )),
        (7, Kind::Mole) => Some(format!(
            " = {}, {}",
            set_units(1.0 * litres * weight / calc_value, "L")?,
            set_units(1.0 * litres * weight, "g")?
        )),
        _ => None,
    }
}

/// A weight is needed to go between mass and moles.
fn needs_weight(input: &MolarityInput) -> bool {
    let (first, second) = (input.unit1, input.unit2);
    (first.is_mass() && second.is_mole())
        || (second.is_mass() && first.is_mole())
        || (first == Unit::PercentWv && second.is_mole())
}

/// A sequence is read for the nucleotide weights, and its length is the weight's length: the sequence as it is kept,
/// and the weight field after it.
pub fn sequence_and_length(input: &MolarityInput) -> (String, String) {
    if !input.weight.takes_sequence() {
        return (String::new(), input.mol_wt.to_string());
    }
    let sequence = clean_sequence(input.sequence, input.weight == Weight::Rna);
    let mol_wt = if sequence.is_empty() {
        input.mol_wt.to_string()
    } else {
        sequence.len().to_string()
    };
    (sequence, mol_wt)
}

/// What the volumes of the first and the second quantity are per, as numbers of litres' worth.
fn volume_factors(input: &MolarityInput) -> (f64, f64) {
    match (input.per1, input.unit1) {
        (Some(per1), _) => (per1.factor(), input.per2.factor()),
        (None, Unit::PercentWv) => (1e1, input.per2.factor()),
        (None, _) => (1.0, 1.0),
    }
}

/// The first quantity in the second unit: through the molecular weight when it goes between mass and moles.
fn converted(input: &MolarityInput, first: f64, weight: f64, volumes: (f64, f64)) -> f64 {
    let (conv1, conv2) = (input.unit1.factor(), input.unit2.factor());
    let (vol_conv1, vol_conv2) = volumes;
    if needs_weight(input) && !input.unit2.is_mass() {
        1.0 * first / weight * conv1 / conv2 * vol_conv1 / vol_conv2
    } else if input.unit2.is_mass() && input.unit1.is_mole() {
        1.0 * first * weight * conv1 / conv2 * vol_conv1 / vol_conv2
    } else {
        1.0 * first * conv1 / conv2 * vol_conv1 / vol_conv2
    }
}

/// Which second line goes with the units (the key: 4 for grams, 2 for moles, 1 for a volume, added up), and the first
/// quantity as grams or moles per litre that it is worked out from (`calc_value`).
fn second_line_key(input: &MolarityInput, first: f64, weight: f64, vol_conv1: f64) -> (u8, f64) {
    let (first_unit, second_unit) = (input.unit1, input.unit2);
    let mut key = 0;
    let mut calc_value = first * first_unit.factor();
    if first_unit.is_mass() || first_unit == Unit::PercentWv || second_unit.is_mass() {
        key += 4;
    }
    if first_unit.is_mole() || second_unit.is_mole() {
        key += 2;
    }
    if first_unit.is_mole() && second_unit.is_mass() {
        calc_value *= weight;
    }
    if input.per1.is_some() || first_unit == Unit::PercentWv {
        key += 1;
        calc_value *= vol_conv1;
    }
    (key, calc_value)
}

/// Convert the first quantity to the second unit, and work out the second line.
pub fn molarity(input: &MolarityInput) -> MolarityOutput {
    let (sequence, mol_wt) = sequence_and_length(input);
    let mut out = MolarityOutput {
        mol_wt: mol_wt.clone(),
        sequence: sequence.clone(),
        needs_weight: needs_weight(input),
        per_volume: input.per1.is_some() || input.unit1 == Unit::PercentWv,
        ..MolarityOutput::default()
    };
    let (Some(first), Some(length), Some(amount)) =
        (number(input.first), number(&mol_wt), number(input.amount))
    else {
        return out;
    };
    if length < 1.0 {
        return out;
    }

    let weight = molecular_weight(input.weight, &sequence, length);
    let volumes = volume_factors(input);
    out.second = g3(converted(input, first, weight, volumes));
    let (key, calc_value) = second_line_key(input, first, weight, volumes.0);
    out.key = key;
    out.detail =
        second_line(key, input.amount_unit, amount, calc_value, weight).unwrap_or_default();
    out
}
