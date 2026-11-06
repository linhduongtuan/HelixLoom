//! The reactions the Molecular Reaction Calculator knows, and the values it starts with. All of it is from AppMain.tcl:
//! the molar ratio reactions (`info(calc_molar_raction_dict)`, L67338, 0x31222b), the master mix (`calc_mm_raction_dict`,
//! L67339) and the first values (`calc_defaults`, L67340).
//!
//! The original also keeps a list of fragments with each molar ratio reaction (their names, sizes, concentrations and
//! ratios), but never loads it while "Preserve Values Between Reactions" is on, which it always is (the check box that
//! would turn it off is commented out), so choosing a reaction only sets the conditions. Those lists are not here.

use super::reaction::DnaUnits;

/// The conditions of a molar ratio reaction.
#[derive(Clone, Copy, Debug)]
pub struct MolarPreset {
    pub name: &'static str,
    /// The volume of the reaction, in ul.
    pub volume: &'static str,
    pub final_amount: &'static str,
    pub units: DnaUnits,
    /// The smallest volume that can be pipetted, in ul.
    pub min_pipette: &'static str,
    /// How concentrated the enzyme mix is (X).
    pub enzyme_x: &'static str,
    pub notes: &'static str,
}

pub const MOLAR_PRESETS: [MolarPreset; 8] = [
    MolarPreset {
        name: "HiFi 2-3 fragment",
        volume: "4",
        final_amount: "40",
        units: DnaUnits::Fmole,
        min_pipette: "0.3",
        enzyme_x: "2",
        notes: "6 to 40 fmole per 4ul. linear:equimolar, insertion:3-to-1,<200 bp ins:5-to-1",
    },
    MolarPreset {
        name: "HiFi 4-6 fragment",
        volume: "4",
        final_amount: "60",
        units: DnaUnits::Fmole,
        min_pipette: "0.3",
        enzyme_x: "2",
        notes: "40 to 100 fmol per 4ul. linear:equimolar, insertion:3-to-1,<200 bp ins:5-to-1",
    },
    MolarPreset {
        name: "Gibson 2-3 fragment",
        volume: "4",
        final_amount: "40",
        units: DnaUnits::Fmole,
        min_pipette: "0.3",
        enzyme_x: "1.333",
        notes: "6 to 40 fmole per 4ul. linear:equimolar, insertion:3-to-1,<200 bp ins:5-to-1",
    },
    MolarPreset {
        name: "Gibson 4-6 fragment",
        volume: "4",
        final_amount: "60",
        units: DnaUnits::Fmole,
        min_pipette: "0.3",
        enzyme_x: "1.333",
        notes: "40 to 200 fmole per 4ul. linear:equimolar, insertion:3-to-1,<200 bp ins:5-to-1",
    },
    MolarPreset {
        name: "SapTrap",
        volume: "2.5",
        final_amount: "30",
        units: DnaUnits::Nm,
        min_pipette: "0.3",
        enzyme_x: "1.25",
        notes: "oligos:plasmid at 3:1",
    },
    MolarPreset {
        name: "Golden Gate",
        volume: "10",
        final_amount: "6",
        units: DnaUnits::Nm,
        min_pipette: "0.3",
        enzyme_x: "5",
        notes: "adjust final conc. to 3 nM each fragment",
    },
    MolarPreset {
        name: "T4 Ligation",
        volume: "10",
        final_amount: "10",
        units: DnaUnits::NgPerUl,
        min_pipette: "0.3",
        enzyme_x: "5",
        notes: "DNA should be 1-10ng/ul. V:I reatio should be between 1:1 and 1:10.",
    },
    MolarPreset {
        name: "NEB Quick Ligation",
        volume: "20",
        final_amount: "10",
        units: DnaUnits::NgPerUl,
        min_pipette: "0.3",
        enzyme_x: "2",
        notes: "DNA should be 1-10ng/ul. V:I reatio should be between 1:2 and 1:6.",
    },
];

/// The choices of the enzyme stock concentration menu, with the X each one sets (`DNA Only` is 1: no enzyme mix).
pub const ENZYME_STOCKS: [(&str, &str); 6] = [
    ("10x (1:9)", "10"),
    ("5x (1:4)", "5"),
    ("2x (1:1)", "2"),
    ("1.33x (3:1)", "1.333"),
    ("1.25x (4:1)", "1.25"),
    ("DNA Only", "1"),
];

/// A master mix: the name and the volume (ul) of each component in one reaction.
pub struct MixPreset {
    pub name: &'static str,
    pub components: &'static [(&'static str, &'static str)],
}

pub const MIX_PRESETS: [MixPreset; 1] = [MixPreset {
    name: "Phusion",
    components: &[
        ("primer", ".25"),
        ("primer", ".25"),
        ("enzyme", ".5"),
        ("dNTP", "1"),
        ("water", "38"),
        ("buffer", "10"),
    ],
}];

/// The values the window opens with (`calc_defaults`).
pub struct Defaults {
    pub molar_reaction: &'static str,
    pub mix_reaction: &'static str,
    pub reactions: &'static str,
    pub extra: &'static str,
    pub molarity_first: &'static str,
    pub molarity_mol_wt: &'static str,
    pub molarity_amount: &'static str,
    pub stock: &'static str,
    pub final_conc: &'static str,
    pub final_volume: &'static str,
}

pub const DEFAULTS: Defaults = Defaults {
    molar_reaction: "HiFi 2-3 fragment",
    mix_reaction: "Phusion",
    reactions: "1",
    extra: "0",
    molarity_first: "1",
    molarity_mol_wt: "1000",
    molarity_amount: "1",
    stock: "1",
    final_conc: "1",
    final_volume: "1",
};

/// The enzyme stock label that matches an X, if any of the menu's choices does.
pub fn enzyme_label(x: &str) -> Option<&'static str> {
    ENZYME_STOCKS
        .iter()
        .find(|(_, value)| *value == x)
        .map(|(label, _)| *label)
}
