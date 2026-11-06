//! The stock dilution: how much of a stock to take, and how much water to add, for a final volume of a final
//! concentration. Re-written from `calc_stock_dilution_update` (AppMain.tcl:35540, 0x19e731) with the arithmetic in the
//! same order.

use super::format::{number, tcl_double};

/// A concentration unit (`X M mM uM nM pM`): X is for what is sold as ten times, five times and so on concentrated.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Conc {
    X,
    M,
    Mm,
    Um,
    Nm,
    Pm,
}

impl Conc {
    pub const ALL: [Conc; 6] = [Conc::X, Conc::M, Conc::Mm, Conc::Um, Conc::Nm, Conc::Pm];

    pub fn label(self) -> &'static str {
        match self {
            Conc::X => "X",
            Conc::M => "M",
            Conc::Mm => "mM",
            Conc::Um => "uM",
            Conc::Nm => "nM",
            Conc::Pm => "pM",
        }
    }

    /// Where it stands in the list `M mM uM nM pM X` that the original counts steps of a thousand in.
    fn place(self) -> i32 {
        match self {
            Conc::M => 0,
            Conc::Mm => 1,
            Conc::Um => 2,
            Conc::Nm => 3,
            Conc::Pm => 4,
            Conc::X => 5,
        }
    }
}

/// The unit of a volume (`L mL uL`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Litres {
    L,
    Ml,
    Ul,
}

impl Litres {
    pub const ALL: [Litres; 3] = [Litres::L, Litres::Ml, Litres::Ul];

    pub fn label(self) -> &'static str {
        match self {
            Litres::L => "L",
            Litres::Ml => "mL",
            Litres::Ul => "uL",
        }
    }
}

/// Which unit was changed last: the original keeps X (a multiple) on both sides, or on neither.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Changed {
    Stock,
    Final,
}

/// The two units after one of them was changed: a stock in X is diluted to a final in X, and a final in X comes from a
/// stock in X; leaving X on one side takes the other side along to the unit that was chosen.
pub fn keep_x_together(stock: Conc, final_conc: Conc, changed: Changed) -> (Conc, Conc) {
    match changed {
        Changed::Stock if stock == Conc::X => (stock, Conc::X),
        Changed::Stock if final_conc == Conc::X => (stock, stock),
        Changed::Final if final_conc == Conc::X => (Conc::X, final_conc),
        Changed::Final if stock == Conc::X => (final_conc, final_conc),
        _ => (stock, final_conc),
    }
}

#[derive(Clone, Copy, Debug)]
pub struct StockInput<'a> {
    pub stock: &'a str,
    pub stock_unit: Conc,
    pub final_conc: &'a str,
    pub final_unit: Conc,
    pub final_volume: &'a str,
    pub volume_unit: Litres,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StockOutput {
    /// How much stock to take (`Error` when the final concentration is above the stock's); empty when not worked out.
    pub stock: String,
    pub water: String,
}

/// A volume in the smallest unit that keeps it at one or more (`0.5 mL` is `500.0 uL`), to hundredths.
fn adjusted(volume: f64, unit: Litres) -> String {
    let (mut value, mut unit) = (volume, unit);
    while value < 1.0 && unit != Litres::Ul {
        value *= 1000.0;
        unit = if unit == Litres::L {
            Litres::Ml
        } else {
            Litres::Ul
        };
    }
    format!(
        "{} {}",
        tcl_double((value * 100.0).round() / 100.0),
        unit.label()
    )
}

/// The three numbers, when all three are there and the stock is not nothing.
fn numbers(input: &StockInput) -> Option<(f64, f64, f64)> {
    let (stock, final_conc, final_volume) = (
        number(input.stock)?,
        number(input.final_conc)?,
        number(input.final_volume)?,
    );
    (stock != 0.0).then_some((stock, final_conc, final_volume))
}

/// What the original does when something was changed: nothing at all while a number is missing, else the units follow
/// X, and then the volumes are worked out. Gives the units that are left and the volumes.
pub fn update(input: &StockInput, changed: Option<Changed>) -> (Conc, Conc, StockOutput) {
    if numbers(input).is_none() {
        return (input.stock_unit, input.final_unit, StockOutput::default());
    }
    let (stock_unit, final_unit) = changed.map_or((input.stock_unit, input.final_unit), |c| {
        keep_x_together(input.stock_unit, input.final_unit, c)
    });
    let input = StockInput {
        stock_unit,
        final_unit,
        ..*input
    };
    (stock_unit, final_unit, stock_dilution(&input))
}

/// How much of the stock to take and how much water to add.
pub fn stock_dilution(input: &StockInput) -> StockOutput {
    let Some((stock, final_conc, final_volume)) = numbers(input) else {
        return StockOutput::default();
    };
    let unit_diff = 10f64.powf(f64::from(
        3 * (input.stock_unit.place() - input.final_unit.place()),
    ));
    let stock_volume = final_volume * final_conc * unit_diff / stock;
    let water_volume = final_volume - stock_volume;
    if water_volume < 0.0 {
        return StockOutput {
            stock: "Error".into(),
            water: String::new(),
        };
    }
    StockOutput {
        stock: adjusted(stock_volume, input.volume_unit),
        water: adjusted(water_volume, input.volume_unit),
    }
}

/// What Copy puts on the clipboard for the dilution (`calc_output`, AppMain.tcl:35616, 0x19f21e, tab 4).
pub fn copy_text(input: &StockInput, out: &StockOutput) -> String {
    format!(
        "[Stock]: {} {}\nFinal: {} {} of {} {}\nAdd : {} of stock\nTo : {} water.\n",
        input.stock,
        input.stock_unit.label(),
        input.final_volume,
        input.volume_unit.label(),
        input.final_conc,
        input.final_unit.label(),
        out.stock,
        out.water
    )
}
