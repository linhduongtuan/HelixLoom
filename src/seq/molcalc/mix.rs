//! The master mix multiplier: the volume of each component of a reaction, multiplied up for a number of reactions with
//! some extra (a share of the whole, or a volume) for what sticks to the pipette. Re-written from `calc_mm_reaction`
//! (AppMain.tcl:35131, 0x199d47) with the arithmetic in the same order.
//!
//! **Deliberate difference:** an extra in ul is divided by the volume of the whole mix (`extra / (single volume x
//! reactions)`), and when the three are typed as whole numbers Tcl divides whole numbers, so the original's extra of
//! 20 ul in 12 reactions of 50 ul comes to 20 / 600 = 0 and nothing is added; it works only when one of them has a decimal
//! point. Here the division is always of real numbers.

use super::format::{number, tcl_double};

/// How many components the table holds.
pub const COMPONENTS: usize = 10;

/// What the person did last: it decides what is worked out again.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MixEdit {
    /// A component's volume was typed: the single reaction volume is the sum of them.
    Component,
    /// The single reaction volume was typed: every component is scaled to it.
    SingleVolume,
    /// The number of reactions, the extra, or the reaction itself.
    Other,
}

/// What the extra is given in.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Extra {
    Percent,
    Ul,
}

impl Extra {
    pub fn label(self) -> &'static str {
        match self {
            Extra::Percent => "%",
            Extra::Ul => "ul",
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct MixInput<'a> {
    /// The volume of each component in one reaction, as typed.
    pub components: &'a [String],
    pub single_volume: &'a str,
    pub reactions: &'a str,
    pub extra: &'a str,
    pub extra_unit: Extra,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MixOutput {
    /// The volumes of the components in one reaction (scaled when the single reaction volume was typed).
    pub components: Vec<String>,
    pub single_volume: String,
    /// The volume of each component for the whole mix; empty for one that cannot be worked out.
    pub volumes: Vec<String>,
    pub total: String,
}

/// The volumes of a master mix.
pub fn master_mix(input: &MixInput, edit: MixEdit) -> MixOutput {
    let n = input.components.len().min(COMPONENTS);
    let mut components: Vec<String> = input.components[..n].to_vec();
    // only what is a number counts
    let sum: f64 = components.iter().filter_map(|c| number(c)).sum();

    let mut single_volume = input.single_volume.to_string();
    if edit == MixEdit::Component {
        single_volume = tcl_double(((100.0 * sum).round() + 0.0) / 100.0);
    }
    let single = number(&single_volume).filter(|s| *s > 0.0);
    let reactions = number(input.reactions).filter(|r| *r > 0.0);
    let extra = number(input.extra);

    let mut volumes = vec![String::new(); n];
    let mut sum2 = 0.0;
    for i in 0..n {
        let (Some(v), Some(single), Some(reactions)) = (number(&components[i]), single, reactions)
        else {
            continue;
        };
        if sum <= 0.0 {
            continue;
        }
        let v = if edit == MixEdit::SingleVolume {
            let scaled = (100.0 * v * single / sum).round() / 100.0;
            components[i] = tcl_double(scaled);
            scaled
        } else {
            v
        };
        let share = match (extra, input.extra_unit) {
            (Some(add), Extra::Percent) => 1.0 + add / 100.0,
            (Some(add), Extra::Ul) => 1.0 + add / (single * reactions),
            (None, _) => 1.0 + 0.0 / 100.0,
        };
        let mixed = (v * reactions * share * 100.0).round() / 100.0;
        volumes[i] = tcl_double(mixed);
        sum2 += mixed;
    }
    MixOutput {
        components,
        single_volume,
        volumes,
        total: tcl_double((100.0 * sum2).round() / 100.0),
    }
}

/// What Copy puts on the clipboard for the master mix (`calc_output`, AppMain.tcl:35616, 0x19f21e, tab 3).
pub fn copy_text(
    reaction: &str,
    names: &[String],
    single_volume: &str,
    reactions: &str,
    out: &MixOutput,
) -> String {
    let mut text = format!(
        "Reaction: {reaction}\nSingle Reaction Volume: {single_volume}\nReaction Number: {reactions}\n"
    );
    for (name, volume) in names.iter().zip(&out.volumes) {
        if !volume.is_empty() {
            text.push_str(&format!("{name}\t{volume}\n"));
        }
    }
    text
}
