//! Floating tool windows and the assembler tab. Each file adds an `impl App` block for one tool.
pub mod assembly;
pub mod calculator;
pub mod contig_view;
pub mod digest;
pub mod domains;
pub mod enzlist;
pub mod enzymes;
pub mod features;
pub mod fetch;
pub mod find;
#[cfg(test)]
mod gff_tests;
pub mod molcalc;
pub mod motifs;
pub mod orfs;
pub mod phylocnn;
pub mod primer3;
#[cfg(test)]
mod primer3_tests;
pub mod programs;
pub mod reads;
pub mod sitetools;
pub mod species;
pub mod textmap;
pub mod tools;
#[cfg(test)]
mod tools_pcr_tests;
pub mod tree;
pub mod treetools;
pub mod update;
#[cfg(test)]
mod update_tests;
