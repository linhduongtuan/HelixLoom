//! Colours for bases, residues and the letters of a trace that stay apart for people who see colour differently.
//!
//! The nucleotides and the trace letters are taken from Okabe and Ito's palette (a bluish green, a sky or deeper blue, a yellow or
//! orange and a vermillion), softened for the cells that hold letters. The eight groups of amino acids need more colours than that
//! palette has, so they were found by a search for the largest smallest distance (CIE Lab) between any two of the eight, as seen with
//! normal vision and with each of three kinds of colour blindness (Machado, Oliveira and Fernandes 2009, full strength), with the colours kept
//! light enough for black letters and weak enough not to be loud. The tests in `tests.rs` repeat that measurement, so a colour
//! changed here that spoils it fails them. The groups are the ones of the standard colouring (hydrophobic, positive, negative, polar,
//! cysteine, glycine, proline, aromatic) but not their hues: it is the hues that cannot be told apart.
//!
//! Colour is never the only sign: the letters are written in the cells as well.

use eframe::egui::Color32;

const fn rgb(r: u8, g: u8, b: u8) -> Color32 {
    Color32::from_rgb(r, g, b)
}

/// The colour of a trace letter (or the letter that a trace calls) for base `b`.
pub fn trace(b: u8, dark: bool) -> Color32 {
    match (b.to_ascii_uppercase(), dark) {
        (b'A', _) => rgb(0, 158, 115),
        (b'C', false) => rgb(0, 114, 178),
        (b'C', true) => rgb(86, 180, 233),
        (b'G', false) => Color32::from_gray(0x22),
        (b'G', true) => Color32::from_gray(0xe0),
        (b'T', _) => rgb(213, 94, 0),
        _ => rgb(204, 121, 167),
    }
}

/// The background of a cell that holds nucleotide `b`, if it has one.
pub fn base(b: u8, dark: bool) -> Option<Color32> {
    Some(match (b.to_ascii_uppercase(), dark) {
        (b'A', false) => rgb(0x4d, 0xbb, 0x9d),
        (b'C', false) => rgb(0x89, 0xca, 0xf0),
        (b'G', false) => rgb(0xf4, 0xec, 0x7b),
        (b'T' | b'U', false) => rgb(0xe2, 0x8e, 0x4d),
        (b'A', true) => rgb(0x00, 0x8e, 0x68),
        (b'C', true) => rgb(0x00, 0x67, 0xa0),
        (b'G', true) => rgb(0xd8, 0xcd, 0x3b),
        (b'T' | b'U', true) => rgb(0xc0, 0x55, 0x00),
        _ => return None,
    })
}

/// The eight groups of amino acids, in the order of their colours below.
const GROUPS: [&[u8]; 8] = [
    b"AILMFWV", // hydrophobic
    b"KR",      // positive
    b"DE",      // negative
    b"NQST",    // polar
    b"C",       // cysteine
    b"G",       // glycine
    b"P",       // proline
    b"HY",      // aromatic
];

const RESIDUES_LIGHT: [Color32; 8] = [
    rgb(0x57, 0xd0, 0xff),
    rgb(0xf9, 0xa7, 0xc3),
    rgb(0xb9, 0x89, 0xc0),
    rgb(0x44, 0xa7, 0x84),
    rgb(0xd5, 0xf0, 0xc6),
    rgb(0x89, 0xaf, 0x6b),
    rgb(0xeb, 0xde, 0x95),
    rgb(0xb7, 0xf0, 0xff),
];

/// The same colours at 65 per cent of their strength, for a dark page.
const RESIDUES_DARK: [Color32; 8] = [
    rgb(0x39, 0x87, 0xa6),
    rgb(0xa2, 0x6d, 0x7f),
    rgb(0x78, 0x59, 0x7d),
    rgb(0x2c, 0x6d, 0x56),
    rgb(0x8a, 0x9c, 0x81),
    rgb(0x59, 0x72, 0x46),
    rgb(0x99, 0x90, 0x61),
    rgb(0x77, 0x9c, 0xa6),
];

/// The background of a cell that holds amino acid `aa`, if it has one.
pub fn residue(aa: u8, dark: bool) -> Option<Color32> {
    let aa = aa.to_ascii_uppercase();
    let group = GROUPS.iter().position(|letters| letters.contains(&aa))?;
    Some(if dark {
        RESIDUES_DARK[group]
    } else {
        RESIDUES_LIGHT[group]
    })
}
