//! Floating tool windows and the assembler tab. Each file adds an `impl App` block for one tool.
pub mod assembly;
pub mod calculator;
pub mod contig_view;
pub mod digest;
pub mod enzymes;
pub mod features;
pub mod find;
pub mod molcalc;
pub mod orfs;
pub mod phylocnn;
pub mod primer3;
pub mod reads;
pub mod tools;
#[cfg(test)]
mod tools_pcr_tests;
pub mod tree;
pub mod treetools;
