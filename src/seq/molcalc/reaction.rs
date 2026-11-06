//! The molar ratio reaction of the Molecular Reaction Calculator: how much of each DNA fragment to pipette so that the
//! fragments meet in the ratio asked for (Gibson, HiFi, Golden Gate and ligation reactions). Re-written from
//! `calc_molar_reaction` (AppMain.tcl:34777, 0x195d0d); the arithmetic is done in the same order, so the numbers are the
//! same to the last digit.

use super::format::{fixed2, number, tcl_double};

/// How many fragments the table holds.
pub const FRAGMENTS: usize = 10;

/// What a kilobase pair of double-stranded DNA weighs, in nanograms per picomole (617.96 g/mol for each base pair), and
/// what the two ends add (36.04): `size * 617.96 + 36.04`.
const NG_PER_PMOLE_PER_KB: f64 = 617.96;
const ENDS: f64 = 36.04;
/// Water that is less than this share of the DNA volume is left out of the smallest pipette volume.
const WATER_SHARE: f64 = 0.05;
/// What the smallest volume starts from: no reaction is as big as this many ul.
const SMALLEST_START: f64 = 10000.0;

/// What the final DNA is given in: the menu of the original (`nM fmole ng ng/ul`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DnaUnits {
    Nm,
    Fmole,
    Ng,
    NgPerUl,
}

impl DnaUnits {
    pub const ALL: [DnaUnits; 4] = [
        DnaUnits::Nm,
        DnaUnits::Fmole,
        DnaUnits::Ng,
        DnaUnits::NgPerUl,
    ];

    pub fn label(self) -> &'static str {
        match self {
            DnaUnits::Nm => "nM",
            DnaUnits::Fmole => "fmole",
            DnaUnits::Ng => "ng",
            DnaUnits::NgPerUl => "ng/ul",
        }
    }

    /// The amount is a mass (ng, or ng per ul) and not moles.
    fn in_ng(self) -> bool {
        matches!(self, DnaUnits::Ng | DnaUnits::NgPerUl)
    }
}

/// One row of the table as typed: any of the fields may be empty or not a number.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Fragment {
    pub name: String,
    pub size_kb: String,
    /// Concentration of the stock, in ng/ul.
    pub conc: String,
    pub ratio: String,
}

/// Everything the calculation reads.
#[derive(Clone, Copy, Debug)]
pub struct MolarInput<'a> {
    pub fragments: &'a [Fragment],
    /// The volume of the whole reaction, in ul.
    pub reaction_volume: &'a str,
    pub final_amount: &'a str,
    pub units: DnaUnits,
    /// The smallest volume that can be pipetted, in ul.
    pub min_pipette: &'a str,
    /// How concentrated the enzyme mix is (10 for 10x, 2 for 2x); 1 or less for DNA only.
    pub enzyme_x: &'a str,
    pub add_water: bool,
}

/// What comes out for one row.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RowOut {
    /// The ratio: as typed, or 1 when the row has a size and a concentration but no ratio.
    pub ratio: String,
    /// The volume to pipette, in ul (two decimals); empty for a row that cannot be used.
    pub volume: String,
    /// How much DNA that is: fmoles, or ng for the units in ng.
    pub amount: String,
}

/// What comes out in all, as the original writes it: empty when it cannot be worked out.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MolarOutput {
    pub rows: Vec<RowOut>,
    /// What the row amounts are in: "fmoles" or "ng".
    pub amount_label: &'static str,
    pub water: String,
    pub intermediate: String,
    /// The DNA mix volume, with how far the DNA is off the volume when it is.
    pub dna_volume: String,
    pub enzyme_volume: String,
    /// The DNA is too dilute for the volume: the original colours the label of the final amount red.
    pub too_dilute: bool,
}

/// What Copy puts on the clipboard (`calc_output`, AppMain.tcl:35616, 0x19f21e, tab 2).
pub fn copy_text(reaction: &str, fragments: &[Fragment], out: &MolarOutput, notes: &str) -> String {
    let mut text =
        format!("Reaction: {reaction}\nFragment Name\tSize (kb)\tng/ul\tratio\tVolume\n");
    for (f, row) in fragments.iter().zip(&out.rows) {
        if !row.volume.is_empty() {
            text.push_str(&format!(
                "{}\t{}\t{}\t{}\t{} ul\n",
                f.name, f.size_kb, f.conc, row.ratio, row.volume
            ));
        }
    }
    text.push_str(&format!(
        "Water: {} ul\n\nDNA mix: {} ul\nEnzyme+Buffer mix: {} ul\nNotes: {notes}",
        out.water, out.dna_volume, out.enzyme_volume
    ));
    text
}

/// A row that can be used, with what it says.
#[derive(Clone, Copy)]
struct Used {
    size: f64,
    conc: f64,
    ratio: f64,
}

/// The table as typed: the rows as they come out, the rows that can be used, and their ratios added up.
struct Table {
    rows: Vec<RowOut>,
    used: Vec<Option<Used>>,
    ratio_sum: f64,
}

/// A row is used when its size, concentration and ratio are numbers that multiply to more than nothing; a ratio left
/// empty is 1 once the row has a size and a concentration.
fn read_table(fragments: &[Fragment]) -> Table {
    let mut rows: Vec<RowOut> = fragments
        .iter()
        .map(|f| RowOut {
            ratio: f.ratio.clone(),
            ..RowOut::default()
        })
        .collect();
    let mut used = vec![None; fragments.len()];
    let mut ratio_sum = 0.0;
    for (i, f) in fragments.iter().enumerate() {
        let (size, conc) = (number(&f.size_kb), number(&f.conc));
        if let (Some(s), Some(c)) = (size, conc) {
            if s * c > 0.0 && f.ratio.trim().is_empty() {
                rows[i].ratio = "1".into();
            }
        }
        if let (Some(size), Some(conc), Some(ratio)) = (size, conc, number(&rows[i].ratio)) {
            if size * conc * ratio > 0.0 {
                ratio_sum += ratio;
                used[i] = Some(Used { size, conc, ratio });
            }
        }
    }
    Table {
        rows,
        used,
        ratio_sum,
    }
}

/// The volume of each fragment's stock (ul) and how much DNA that is, for the rows that are used (0 for the others).
struct Volumes {
    vol: Vec<f64>,
    /// fmoles, or ng for the units in ng.
    amt: Vec<f64>,
    /// All the DNA, in ul.
    dna_sum: f64,
}

impl Volumes {
    /// Every row that is used, times `factor`.
    fn scale(&mut self, used: &[Option<Used>], factor: f64) {
        for (i, u) in used.iter().enumerate() {
            if u.is_some() {
                self.vol[i] *= factor;
                self.amt[i] *= factor;
            }
        }
    }

    /// The smallest volume of a row that is used.
    fn smallest(&self, used: &[Option<Used>]) -> f64 {
        let mut smallest = SMALLEST_START;
        for (i, u) in used.iter().enumerate() {
            if u.is_some() {
                smallest = smallest.min(self.vol[i]);
            }
        }
        smallest
    }
}

/// The volume of a fragment's stock and, for nM and fmoles, how much DNA that is. For the units in ng it is only the
/// share of the volume (the ratio), which `fragment_volumes` scales to the ng asked for.
fn row_volume(units: DnaUnits, u: &Used, ratio_sum: f64, volume: f64, amount: f64) -> (f64, f64) {
    let grams = u.size * NG_PER_PMOLE_PER_KB + ENDS;
    match units {
        DnaUnits::Nm => {
            let vol = 1.0 * amount / 1000.0 * volume * u.ratio / ratio_sum * grams / u.conc;
            (vol, 1.0 * 1000.0 * vol * u.conc / grams)
        }
        DnaUnits::Fmole => {
            let vol = 1.0 * amount / 1000.0 * u.ratio / ratio_sum * grams / u.conc;
            (vol, 1.0 * 1000.0 * vol * u.conc / grams)
        }
        DnaUnits::Ng | DnaUnits::NgPerUl => (1.0 * u.ratio / ratio_sum * u.size / u.conc, 0.0),
    }
}

/// The volumes of all the fragments that are used.
fn fragment_volumes(input: &MolarInput, table: &Table, volume: f64, amount: f64) -> Volumes {
    let n = table.used.len();
    let mut v = Volumes {
        vol: vec![0.0; n],
        amt: vec![0.0; n],
        dna_sum: 0.0,
    };
    let mut ng_sum = 0.0;
    for (i, u) in table.used.iter().enumerate() {
        let Some(u) = u else {
            continue;
        };
        (v.vol[i], v.amt[i]) = row_volume(input.units, u, table.ratio_sum, volume, amount);
        ng_sum += v.vol[i] * u.conc;
        v.dna_sum += v.vol[i];
    }
    if input.units.in_ng() {
        // the volumes above only give the ratio: scale them to the total ng asked for
        v.dna_sum = 0.0;
        for (i, u) in table.used.iter().enumerate() {
            let Some(u) = u else {
                continue;
            };
            v.vol[i] = if input.units == DnaUnits::Ng {
                amount / ng_sum * v.vol[i]
            } else {
                volume * amount / ng_sum * v.vol[i]
            };
            v.amt[i] = u.conc * v.vol[i];
            v.dna_sum += v.vol[i];
        }
    }
    v
}

/// The enzyme mix takes its share of the reaction (1/X of it), and the DNA mix the rest: the volume of the DNA mix as a
/// number and as it is written, and the volume of the enzyme mix as it is written.
fn enzyme_split(input: &MolarInput, volume: f64) -> (f64, String, String) {
    match number(input.enzyme_x).filter(|x| *x > 1.0) {
        Some(x) => {
            let enzyme = ((100.0 * volume) / x).round() / 100.0;
            let dna = ((100.0 * volume) * (1.0 - 1.0 / x)).round() / 100.0;
            (dna, tcl_double(dna), tcl_double(enzyme))
        }
        None => (volume, input.reaction_volume.to_string(), "0".to_string()),
    }
}

/// The DNA mix: its volume, the water that makes it up, and how much of the volume the DNA is (a percentage, 100 when
/// it fits).
struct Mix {
    dna_inter: f64,
    water: f64,
    percent: i64,
}

impl Mix {
    /// Water makes up the DNA mix volume; with too little DNA, or without water, the volumes are scaled to fill it.
    fn filled(add_water: bool, dna_inter: f64, v: &mut Volumes, used: &[Option<Used>]) -> Mix {
        if add_water && dna_inter >= v.dna_sum {
            return Mix {
                dna_inter,
                water: dna_inter - v.dna_sum,
                percent: 100,
            };
        }
        let percent = (100.0 * dna_inter / v.dna_sum) as i64;
        let scale = dna_inter / v.dna_sum;
        v.dna_sum = dna_inter;
        v.scale(used, scale);
        Mix {
            dna_inter,
            water: 0.0,
            percent,
        }
    }

    /// The smallest volume to pipette: the water counts when it is at least 5 % of the DNA, else it is left out.
    fn smallest(&mut self, v: &Volumes, used: &[Option<Used>]) -> f64 {
        let mut smallest = v.smallest(used);
        if self.water < smallest && self.percent == 100 {
            if self.water >= WATER_SHARE * v.dna_sum {
                smallest = smallest.min(self.water);
            } else {
                self.water = 0.0;
                self.percent = (100.0 * self.dna_inter / v.dna_sum) as i64;
            }
        }
        smallest
    }

    /// The DNA mix volume as it is written (`volume`), with how far the DNA is off it when it is.
    fn dna_text(&self, volume: String) -> String {
        if self.percent < 100 {
            format!(
                "{volume} -Input DNA concentrations are too low ({} %)",
                self.percent
            )
        } else if self.percent > 100 {
            format!("{volume} ({} %)", self.percent)
        } else {
            volume
        }
    }
}

/// Nothing is worked out: only the ratios that were filled in.
fn only_rows(rows: Vec<RowOut>) -> MolarOutput {
    MolarOutput {
        rows,
        ..MolarOutput::default()
    }
}

/// The volumes of the fragments (and of the water and the enzyme mix) of a molar ratio reaction.
pub fn molar_reaction(input: &MolarInput) -> MolarOutput {
    let mut table = read_table(&input.fragments[..input.fragments.len().min(FRAGMENTS)]);
    let (Some(volume), Some(amount)) = (number(input.reaction_volume), number(input.final_amount))
    else {
        return only_rows(table.rows);
    };
    if table.ratio_sum == 0.0 || volume * amount <= 0.0 || amount < 0.0 {
        return only_rows(table.rows);
    }

    let (dna_inter, dna_text, enzyme_text) = enzyme_split(input, volume);
    let mut v = fragment_volumes(input, &table, volume, amount);
    let mut mix = Mix::filled(input.add_water, dna_inter, &mut v, &table.used);
    let smallest = mix.smallest(&v, &table.used);
    if let Some(min_pipette) = number(input.min_pipette) {
        if smallest < min_pipette && smallest > 0.0 {
            // everything is made bigger until the smallest volume can be pipetted
            let scale = min_pipette / smallest;
            v.scale(&table.used, scale);
            mix.dna_inter *= scale;
            mix.water *= scale;
        }
    }

    for (i, u) in table.used.iter().enumerate() {
        if u.is_some() {
            table.rows[i].volume = fixed2(v.vol[i]);
            table.rows[i].amount = fixed2(v.amt[i]);
        }
    }
    MolarOutput {
        rows: table.rows,
        amount_label: if input.units.in_ng() { "ng" } else { "fmoles" },
        water: fixed2(mix.water),
        intermediate: fixed2(mix.dna_inter),
        dna_volume: mix.dna_text(dna_text),
        enzyme_volume: enzyme_text,
        too_dilute: mix.percent < 100,
    }
}
