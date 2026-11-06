//! The Molecular Reaction Calculator of the original (Tools ▸ Molecular Reaction Calculator, `mol_calc`, AppMain.tcl:34255,
//! 0x18dde0): the volumes of a molar ratio reaction (Gibson, HiFi, Golden Gate, ligation), a master mix multiplier,
//! conversions between mass, moles and concentration, and a stock dilution. Each part was re-written from the original's own
//! procedure and checked against it: `golden/*.txt` hold what the original's procedures printed (run unchanged with the
//! system `tclsh`, 8.5.9) for random inputs, and the tests here give the same inputs to the Rust code.
//!
//! Not built: the lookup of a compound's molecular weight by name (it asks PubChem over the internet), the recipe
//! multiplier (it reads the clipboard and asks PubChem), saving a reaction of your own, dropping gel bands into the
//! table, and the QR code of the result.

pub mod expr;
pub mod format;
pub mod mix;
pub mod molarity;
pub mod presets;
pub mod reaction;
pub mod stock;

#[cfg(test)]
mod tests;
