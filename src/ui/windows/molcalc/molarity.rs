//! The Molarity tab: conversions between mass, moles, concentration and volume (`calc_molarity`) and, below them, the
//! stock dilution (`calc_stock_dilution_update`).

use super::{combo, number_field, result};
use crate::seq::molcalc::molarity::{
    molarity, sequence_and_length, Amount, MolarityInput, MolarityOutput, Per, Unit, Weight,
};
use crate::seq::molcalc::presets::DEFAULTS;
use crate::seq::molcalc::stock::{
    copy_text, stock_dilution, update, Changed, Conc, Litres, StockInput,
};
use eframe::egui::{self, RichText, Ui};

pub struct MolarityTab {
    first: String,
    unit1: Unit,
    per1: Option<Per>,
    unit2: Unit,
    per2: Per,
    weight: Weight,
    mol_wt: String,
    sequence: String,
    amount: String,
    amount_unit: Amount,
    stock: String,
    stock_unit: Conc,
    final_conc: String,
    final_unit: Conc,
    final_volume: String,
    volume_unit: Litres,
}

impl Default for MolarityTab {
    /// The first values of the original (`calc_defaults`).
    fn default() -> MolarityTab {
        MolarityTab {
            first: DEFAULTS.molarity_first.into(),
            unit1: Unit::Mmole,
            per1: Some(Per::L),
            unit2: Unit::G,
            per2: Per::L,
            weight: Weight::GramsPerMole,
            mol_wt: DEFAULTS.molarity_mol_wt.into(),
            sequence: String::new(),
            amount: DEFAULTS.molarity_amount.into(),
            amount_unit: Amount::Ml,
            stock: DEFAULTS.stock.into(),
            stock_unit: Conc::Mm,
            final_conc: DEFAULTS.final_conc.into(),
            final_unit: Conc::Mm,
            final_volume: DEFAULTS.final_volume.into(),
            volume_unit: Litres::L,
        }
    }
}

/// The name of the box for a sequence, so that the keyboard can be given to it.
const SEQUENCE_BOX: &str = "molarity_sequence";

/// The volume a quantity is per, or none.
fn per_label(per: Option<Per>) -> &'static str {
    per.map_or("-", Per::label)
}

/// A row of the dilution: a label, a number and its unit of concentration; true when the unit was chosen anew.
fn concentration_row(
    ui: &mut Ui,
    label: &str,
    id: &str,
    (value, unit): (&mut String, &mut Conc),
) -> bool {
    ui.label(label);
    number_field(ui, value, 64.0);
    let chosen = combo(ui, id, unit, &Conc::ALL, Conc::label);
    ui.end_row();
    chosen
}

impl MolarityTab {
    fn input(&self) -> MolarityInput<'_> {
        MolarityInput {
            first: &self.first,
            unit1: self.unit1,
            per1: self.per1,
            unit2: self.unit2,
            per2: self.per2,
            weight: self.weight,
            mol_wt: &self.mol_wt,
            sequence: &self.sequence,
            amount: &self.amount,
            amount_unit: self.amount_unit,
        }
    }

    fn dilution_input(&self) -> StockInput<'_> {
        StockInput {
            stock: &self.stock,
            stock_unit: self.stock_unit,
            final_conc: &self.final_conc,
            final_unit: self.final_unit,
            final_volume: &self.final_volume,
            volume_unit: self.volume_unit,
        }
    }

    pub fn show(&mut self, ui: &mut Ui) {
        let out = molarity(&self.input());
        self.conversion_row(ui, &out);
        if out.needs_weight {
            self.weight_row(ui);
        }
        if out.per_volume {
            self.amount_row(ui, &out);
        }
        // a sequence that is typed is kept as the original keeps it, and its length is the weight's length; read after
        // the boxes were drawn, so that what was typed in this frame is what is kept
        if self.weight.takes_sequence() {
            (self.sequence, self.mol_wt) = sequence_and_length(&self.input());
        }
        ui.separator();
        self.dilution(ui);
    }

    /// `first unit /volume = second unit /volume`.
    fn conversion_row(&mut self, ui: &mut Ui, out: &MolarityOutput) {
        ui.horizontal(|ui| {
            number_field(ui, &mut self.first, 64.0);
            if combo(ui, "mol_unit1", &mut self.unit1, &Unit::FIRST, Unit::label)
                && self.unit1 == Unit::PercentWv
            {
                // a percentage is already per volume
                self.per1 = None;
            }
            let options = [None, Some(Per::L), Some(Per::Ml), Some(Per::Ul)];
            combo(ui, "mol_per1", &mut self.per1, &options, per_label);
            ui.label("=");
            result(ui, &out.second);
            combo(ui, "mol_unit2", &mut self.unit2, &Unit::SECOND, Unit::label);
            if out.per_volume {
                combo(ui, "mol_per2", &mut self.per2, &Per::ALL, Per::label);
            }
        });
    }

    /// The molecular weight: typed in g/mol, or a length or a sequence of DNA or RNA.
    fn weight_row(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            let typed = !self.weight.takes_sequence() || self.sequence.is_empty();
            ui.add_enabled_ui(typed, |ui| number_field(ui, &mut self.mol_wt, 70.0));
            if combo(
                ui,
                "mol_weight",
                &mut self.weight,
                &Weight::ALL,
                Weight::label,
            ) {
                // the original clears the weight for g/mol and starts a length at 1
                self.mol_wt = if self.weight == Weight::GramsPerMole {
                    String::new()
                } else {
                    "1".into()
                };
                self.sequence.clear();
            }
            if self.weight.takes_sequence() {
                ui.add(
                    egui::TextEdit::singleline(&mut self.sequence)
                        .id(egui::Id::new(SEQUENCE_BOX))
                        .hint_text("or the sequence")
                        .desired_width(200.0),
                );
            }
        });
    }

    /// What an amount of the substance comes to: how much of it is in a volume, or what volume holds an amount.
    fn amount_row(&mut self, ui: &mut Ui, out: &MolarityOutput) {
        ui.horizontal(|ui| {
            number_field(ui, &mut self.amount, 64.0);
            combo(
                ui,
                "mol_amount",
                &mut self.amount_unit,
                &Amount::ALL,
                Amount::label,
            );
            ui.label(&out.detail);
        });
    }

    fn dilution(&mut self, ui: &mut Ui) {
        ui.label(RichText::new("Stock dilution").weak());
        let shown = stock_dilution(&self.dilution_input());
        let mut changed = None;
        egui::Grid::new("stock_dilution")
            .num_columns(3)
            .show(ui, |ui| {
                let stock = (&mut self.stock, &mut self.stock_unit);
                if concentration_row(ui, "Stock Concentration", "stock_unit", stock) {
                    changed = Some(Changed::Stock);
                }
                let wanted = (&mut self.final_conc, &mut self.final_unit);
                if concentration_row(ui, "Final Concentration", "final_unit", wanted) {
                    changed = Some(Changed::Final);
                }
                ui.label("Final Volume");
                number_field(ui, &mut self.final_volume, 64.0);
                combo(
                    ui,
                    "volume_unit",
                    &mut self.volume_unit,
                    &Litres::ALL,
                    Litres::label,
                );
                ui.end_row();
                ui.label("Stock");
                result(ui, &shown.stock);
                ui.end_row();
                ui.label("Water");
                result(ui, &shown.water);
                ui.end_row();
            });
        // a stock in X is diluted to a final in X: the other side follows
        let (stock_unit, final_unit, _) = update(&self.dilution_input(), changed);
        (self.stock_unit, self.final_unit) = (stock_unit, final_unit);
        if ui.button("Copy").clicked() {
            let text = copy_text(&self.dilution_input(), &shown);
            ui.ctx().copy_text(text);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ui::headless::frame_input;
    use eframe::egui::{Context, Event, Id};

    /// One frame of the tab, with `typed` typed into the sequence box when it is given the keyboard.
    fn frame(ctx: &Context, tab: &mut MolarityTab, n: u32, typed: Option<&str>) {
        let mut input = frame_input(n);
        if let Some(text) = typed {
            ctx.memory_mut(|m| m.request_focus(Id::new(SEQUENCE_BOX)));
            input.events.push(Event::Text(text.into()));
        }
        let mut output = ctx.run_ui(input, |ui| tab.show(ui));
        output.textures_delta.clear();
    }

    fn dna_tab() -> MolarityTab {
        MolarityTab {
            weight: Weight::SsDna,
            mol_wt: "1".into(),
            ..MolarityTab::default()
        }
    }

    #[test]
    fn a_sequence_typed_into_its_box_is_kept_cleaned_and_its_length_is_the_weights() {
        let ctx = Context::default();
        let mut tab = dna_tab();
        frame(&ctx, &mut tab, 0, None);
        frame(&ctx, &mut tab, 1, Some("acgtx"));
        assert_eq!(tab.sequence, "acgt");
        assert_eq!(tab.mol_wt, "4");
        // more typed goes after it, and the length follows
        frame(&ctx, &mut tab, 2, Some("UA"));
        assert_eq!(tab.sequence, "acgtTA");
        assert_eq!(tab.mol_wt, "6");
    }

    #[test]
    fn the_weight_of_a_kind_of_nucleotide_starts_at_a_length_of_one() {
        let ctx = Context::default();
        let mut tab = MolarityTab::default();
        frame(&ctx, &mut tab, 0, None);
        assert_eq!(
            (tab.weight, tab.mol_wt.as_str()),
            (Weight::GramsPerMole, "1000")
        );
        tab.weight = Weight::Rna;
        tab.mol_wt = "1".into();
        frame(&ctx, &mut tab, 1, None);
        assert_eq!(tab.mol_wt, "1");
        assert!(tab.sequence.is_empty());
    }
}
