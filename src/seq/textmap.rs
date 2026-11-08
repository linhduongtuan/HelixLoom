//! Text Map (Enzymes menu): the sequence as a printable page of text, row by row, with the restriction sites written above the
//! DNA, the second strand, translations in one, three or six frames, an index line, the features and the genes. Built from the
//! original's own text: `enz_text_map` (AppMain.tcl L19401, the dialog) and `do_text_map` (L19730 @0xdec9f, the page), with
//! `linear_enz_sort` (L28745), `translate` (L14518), `get_tags_exact` (L12530), `get_cds` (L14593), `copy_tags` (L12999) and
//! `uppercase_cdna` (L20180). The defaults are the original's (L66936-66955). Where the original stops with a Tcl error the
//! error is kept (see `check` and `Error`), and where it is odd the oddity is kept too: MECHANICS.md section 8, "Text Map".

mod features;
mod genes;
mod highlight;
mod out;
mod tcl;

#[cfg(test)]
mod tests;

pub use features::MapFeature;
use features::{CdsStore, SeqTag};
pub use genes::STANDARD_CODE;
use out::Out;
pub use out::TagInfo;
use std::collections::BTreeMap;
use tcl::{chars, divide, modulo, range, replace, spaces, text};

/// How many characters the left margin of every line is: the line numbers go there.
pub const MARGIN: usize = 7;

/// The kinds of lines the page is made of, in the order the original lists them (`info(enz_text_dialog_order)`, L66936); the
/// person can put them in another order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Frame {
    Index,
    Enzymes,
    Dna,
    Translation,
    SecondStrand,
    Genes,
    Features,
    Space,
}

impl Frame {
    /// The start: index line, enzymes, DNA, translation, second strand, genes, features, space.
    pub const ALL: [Frame; 8] = [
        Frame::Index,
        Frame::Enzymes,
        Frame::Dna,
        Frame::Translation,
        Frame::SecondStrand,
        Frame::Genes,
        Frame::Features,
        Frame::Space,
    ];

    /// What the dialog calls the frame.
    pub fn label(self) -> &'static str {
        match self {
            Frame::Index => "Index Line",
            Frame::Enzymes => "Enzymes",
            Frame::Dna => "DNA",
            Frame::Translation => "Translation",
            Frame::SecondStrand => "2nd Strand",
            Frame::Genes => "Uppercase as Genes",
            Frame::Features => "Features",
            Frame::Space => "Space",
        }
    }
}

/// Where the DNA lines have their numbers (`info(enz_text_dna_numbers)`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LineNumbers {
    None,
    Left,
    Right,
    Both,
}

impl LineNumbers {
    pub const ALL: [LineNumbers; 4] = [
        LineNumbers::None,
        LineNumbers::Left,
        LineNumbers::Right,
        LineNumbers::Both,
    ];

    pub fn label(self) -> &'static str {
        match self {
            LineNumbers::None => "None",
            LineNumbers::Left => "Left",
            LineNumbers::Right => "Right",
            LineNumbers::Both => "Both",
        }
    }

    fn left(self) -> bool {
        matches!(self, LineNumbers::Left | LineNumbers::Both)
    }

    fn right(self) -> bool {
        matches!(self, LineNumbers::Right | LineNumbers::Both)
    }
}

/// How the translation lines read (`info(enz_text_trans)`: "1_1" to "6_3").
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Reading {
    /// 1, 3 or 6 frames.
    pub frames: u8,
    /// Whether an amino acid is written with three letters.
    pub three_letters: bool,
}

impl Reading {
    pub const ALL: [Reading; 6] = [
        Reading::new(1, false),
        Reading::new(3, false),
        Reading::new(6, false),
        Reading::new(1, true),
        Reading::new(3, true),
        Reading::new(6, true),
    ];

    pub const fn new(frames: u8, three_letters: bool) -> Reading {
        Reading {
            frames,
            three_letters,
        }
    }

    /// What the dialog calls it: "3 frame 1 letter".
    pub fn label(self) -> String {
        format!(
            "{} frame {} letter",
            self.frames,
            if self.three_letters { 3 } else { 1 }
        )
    }
}

/// Everything the dialog holds. The start is the original's (L66936-66955).
#[derive(Clone, Debug, PartialEq)]
pub struct Settings {
    pub order: Vec<Frame>,
    /// Bases to a line (`enz_text_width`, 10 to 5000 in the dialog).
    pub width: usize,
    pub show_enzymes: bool,
    /// Blank lines under each group of lines (`enz_text_space`, 1 to 7).
    pub space: usize,
    pub show_space: bool,
    pub numbers: LineNumbers,
    /// Whether the DNA lines carry the colours of the features (`enz_text_copy_highlight`).
    pub copy_highlight: bool,
    pub second_strand: bool,
    pub show_translation: bool,
    pub reading: Reading,
    pub genes: bool,
    pub show_index: bool,
    pub index_char1: char,
    pub index_spacing1: usize,
    pub index_numbers: bool,
    pub index_char2: char,
    pub index_spacing2: usize,
    pub show_features: bool,
    pub show_hidden_features: bool,
    pub crop_features: bool,
}

impl Default for Settings {
    fn default() -> Settings {
        Settings {
            order: Frame::ALL.to_vec(),
            width: 100,
            show_enzymes: false,
            space: 2,
            show_space: true,
            numbers: LineNumbers::Both,
            copy_highlight: true,
            second_strand: false,
            show_translation: false,
            reading: Reading::new(1, false),
            genes: false,
            show_index: true,
            index_char1: '*',
            index_spacing1: 10,
            index_numbers: true,
            index_char2: '*',
            index_spacing2: 5,
            show_features: true,
            show_hidden_features: false,
            crop_features: false,
        }
    }
}

/// The characters the dialog offers for the marks of the index line.
pub const INDEX_CHARS: [char; 8] = [' ', '+', '*', ':', '|', '.', '^', '-'];

/// An enzyme and the places of its sites as `findenzymes` numbers them (from 1, in the order found).
#[derive(Clone, Copy, Debug)]
pub struct EnzymeSites<'a> {
    pub name: &'a str,
    pub sites: &'a [usize],
}

/// The part of the sequence the page is about.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Region {
    /// The whole sequence.
    All,
    /// The bases from `.0` up to but not including `.1`.
    Part(usize, usize),
    /// A selection that runs through the origin of a circle, which the original cannot map: it stops with a Tcl error.
    /// (The selections of this program never do; the page can still be asked for one.)
    #[allow(dead_code)]
    ThroughOrigin,
}

pub struct Request<'a> {
    pub seq: &'a [u8],
    pub circular: bool,
    /// The name of the sequence's window, which heads the page.
    pub title: &'a str,
    /// The line at the top (the original writes the date and the time).
    pub date: &'a str,
    pub region: Region,
    /// Where the selection of the sequence is, whatever part is mapped: the original's text widget has it as a tag, and
    /// where it begins and ends decides how the colours on the DNA lines are put on.
    pub selection: Option<(usize, usize)>,
    /// How many bases a row of the sequence window holds (0 when it is not known). Two oddities of the original depend on it:
    /// the colours on the DNA lines lose their last piece when the sequence is a whole number of rows long, and a feature that
    /// begins at the second base of a row does not show on a line that ends at its twentieth (and so on: see `same_number`).
    pub row_width: usize,
    pub enzymes: &'a [EnzymeSites<'a>],
    pub features: &'a [MapFeature],
    /// The genetic code, 64 letters, in the order of t, c, a, g.
    pub code: &'a str,
    pub settings: &'a Settings,
}

/// What was asked for cannot be made: the Tcl error the original stops with.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Error(pub String);

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

/// A line of DNA: where it is in the text, and the number of its first base as the line shows it (from 1, wrapping at the
/// length of the sequence), and how many bases it has.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DnaRow {
    pub line: usize,
    pub number: usize,
    pub bases: usize,
}

pub struct TextMap {
    pub text: String,
    /// Every tag of the text, the one drawn over the others last. The tags of the features are named `f<number>#` after the
    /// feature's place in the list, the tags of the DNA lines `dnatag` and `ls<number>_1`.
    pub tags: Vec<TagInfo>,
    pub dna_rows: Vec<DnaRow>,
}

/// The part of the dialog that can be wrong before anything is made: the original divides by these.
fn check(s: &Settings) -> Result<(), Error> {
    if s.width == 0 {
        return Err(Error("The line width must be at least 1.".into()));
    }
    if s.show_index && (s.index_spacing1 == 0 || s.index_spacing2 == 0) {
        return Err(Error("divide by zero".into()));
    }
    Ok(())
}

/// The page, or nothing when there is nothing to map (the original returns without a window).
pub fn text_map(req: &Request) -> Result<Option<TextMap>, Error> {
    check(req.settings)?;
    let seq_length = req.seq.len();
    let (start, end) = match req.region {
        Region::All => (0, seq_length),
        Region::Part(a, b) => (a.min(seq_length), b.min(seq_length)),
        Region::ThroughOrigin => {
            return Err(Error("can't read \"seq_len\": no such variable".into()))
        }
    };
    if end <= start {
        return Ok(None);
    }
    let mut run = Run::new(req, start, end);
    run.header();
    run.rows()?;
    Ok(Some(run.finish()))
}

/// What the page is made from while it is made.
struct Run<'a> {
    req: &'a Request<'a>,
    s: &'a Settings,
    /// The bases looked at.
    text: &'a [u8],
    start: usize,
    seq_length: usize,
    /// The translation lines for the whole text, under one another.
    translations: Vec<Vec<char>>,
    genes_line: Vec<char>,
    tags: Vec<SeqTag>,
    cds: BTreeMap<usize, CdsStore>,
    enz_above: bool,
    index_above: bool,
    /// The DNA line before this one: `textline` of the original, which the second strand is made from.
    textline: Option<String>,
    out: Out,
    dna_rows: Vec<DnaRow>,
}

impl<'a> Run<'a> {
    fn new(req: &'a Request<'a>, start: usize, end: usize) -> Run<'a> {
        let s = req.settings;
        let text = &req.seq[start..end];
        let code = req.code.as_bytes();
        let tags = features::sequence_tags(req.features, req.seq.len(), req.circular);
        let position = |f: Frame| s.order.iter().position(|&o| o == f);
        let mut run = Run {
            req,
            s,
            text,
            start,
            seq_length: req.seq.len(),
            translations: Vec::new(),
            genes_line: Vec::new(),
            tags,
            cds: BTreeMap::new(),
            enz_above: position(Frame::Enzymes) < position(Frame::Dna),
            index_above: position(Frame::Index) < position(Frame::Dna),
            textline: None,
            out: Out::new(),
            dna_rows: Vec::new(),
        };
        if s.show_translation {
            run.translations = translations(text, code, s.reading);
        }
        if s.show_features {
            run.cds = features::cds_stores(
                &run.tags,
                req.features,
                req.seq,
                (start, end),
                s.show_hidden_features,
                code,
                run.row(),
            );
        }
        if s.genes {
            run.genes_line = genes::genes_line(text, code);
        }
        run
    }

    /// The date and the name of the sequence, cut into lines that fit.
    fn header(&mut self) {
        self.out.insert(&format!("{}\n", self.req.date), &[]);
        let title = chars(self.req.title);
        let step = self.s.width + 8;
        for u in (0..title.len()).step_by(step) {
            let part = range(&title, u as i64, (u + step - 1) as i64);
            self.out.insert(&format!("{}\n", text(&part)), &[]);
        }
    }

    fn rows(&mut self) -> Result<(), Error> {
        let width = self.s.width;
        let any_shown = self.s.show_enzymes
            || self.s.second_strand
            || self.s.show_translation
            || self.s.show_index
            || self.s.genes
            || self.s.show_features;
        for x in (0..self.text.len()).step_by(width) {
            let linelength = (self.text.len() - x).min(width);
            for &frame in &self.s.order {
                self.frame(frame, x, linelength)?;
            }
            if any_shown {
                self.out.insert("\n", &[]);
            }
        }
        Ok(())
    }

    fn frame(&mut self, frame: Frame, x: usize, linelength: usize) -> Result<(), Error> {
        match frame {
            Frame::Enzymes if self.s.show_enzymes => self.enzymes(x, linelength),
            Frame::Space if self.s.show_space => {
                self.out.insert(&"\n".repeat(self.s.space), &[]);
            }
            Frame::Dna => self.dna(x, linelength),
            Frame::SecondStrand if self.s.second_strand => return self.second_strand(),
            Frame::Translation if self.s.show_translation => self.translation(x, linelength),
            Frame::Index if self.s.show_index => self.index(x, linelength),
            Frame::Genes if self.s.genes => self.genes(x, linelength),
            Frame::Features if self.s.show_features => self.features(x, linelength),
            _ => {}
        }
        Ok(())
    }

    /// How many bases a row of the sequence window has: what the window says, but no more than the sequence (a window with
    /// fewer bases than a row has one row as long as they are). Not known is a row that holds the whole sequence and more.
    fn row(&self) -> usize {
        match self.req.row_width {
            0 => self.seq_length + 1,
            width => width.min(self.seq_length),
        }
    }

    /// Whether the sequence ends where a row of the sequence window ends.
    fn last_row_is_whole(&self) -> bool {
        self.seq_length.is_multiple_of(self.row())
    }

    fn margin(&self) -> String {
        " ".repeat(MARGIN)
    }

    /// The row's part of a line that spans the whole text: `string range line x x+linelength-1`.
    fn part(&self, line: &[char], x: usize, linelength: usize) -> String {
        text(&range(line, x as i64, (x + linelength) as i64 - 1))
    }

    fn dna(&mut self, x: usize, linelength: usize) {
        let textline = String::from_utf8_lossy(&self.text[x..x + linelength]).into_owned();
        let number = seq_len_mod(x + self.start + 1, self.seq_length);
        let numbers = self.s.numbers;
        if numbers.left() {
            let pad = " ".repeat(6usize.saturating_sub(number.to_string().len()));
            self.out.insert(&format!("{pad}{number} "), &[]);
        } else {
            let margin = self.margin();
            self.out.insert(&margin, &[]);
        }
        let line = self.out.line_start().0;
        let ls = format!("ls{number}_1");
        self.out.insert(&textline, &["dnatag", ls.as_str()]);
        self.dna_rows.push(DnaRow {
            line,
            number,
            bases: linelength,
        });
        if self.s.copy_highlight {
            let first = x + self.start;
            let marks = highlight::copy_tags(
                &self.tags,
                self.req.selection,
                (first, first + linelength - 1),
                &textline,
                self.req.features,
                self.seq_length,
                self.last_row_is_whole(),
            );
            for (name, background, from, to) in marks {
                self.out.set_background(&name, background);
                self.out
                    .tag_add(&name, (line, MARGIN + from), (line, MARGIN + to));
            }
        }
        let tail = if numbers.right() {
            let last = seq_len_mod(x + self.start + linelength, self.seq_length);
            format!(" {last}\n")
        } else {
            " \n".to_string()
        };
        // `insert end " 123\n" {}`: with a list of tags of its own, an empty one
        self.out.insert_untagged(&tail);
        self.textline = Some(textline);
    }

    fn second_strand(&mut self) -> Result<(), Error> {
        let Some(line) = &self.textline else {
            return Err(Error("can't read \"textline\": no such variable".into()));
        };
        let other: String = line.chars().map(complement).collect();
        let row = format!("{}{other}\n", self.margin());
        self.out.insert(&row, &[]);
        Ok(())
    }

    fn translation(&mut self, x: usize, linelength: usize) {
        for line in &self.translations {
            let row = format!("{}{}\n", " ".repeat(MARGIN), self.part(line, x, linelength));
            self.out.insert(&row, &[]);
        }
    }

    fn genes(&mut self, x: usize, linelength: usize) {
        let row = format!(
            "{}{}\n",
            " ".repeat(MARGIN),
            self.part(&self.genes_line, x, linelength)
        );
        self.out.insert(&row, &[]);
    }

    /// The restriction sites in this row: a tick under or over the DNA for each, and its place and name written beside it, on
    /// as many lines as it takes to keep the names apart.
    fn enzymes(&mut self, x: usize, linelength: usize) {
        let first = self.start + x + 1;
        let sorted = linear_enz_sort(
            self.req.enzymes,
            first,
            first + linelength - 1,
            self.seq_length,
        );
        let mut tickline = spaces(linelength as i64);
        let mut lines: Vec<Vec<char>> = Vec::new();
        for (site, name) in sorted {
            let column = site as i64 - self.start as i64 - x as i64 - 1;
            let mut i = 0;
            while i < lines.len() && lines[i].len() as i64 >= column {
                i += 1;
            }
            if i == lines.len() {
                lines.push(Vec::new());
            }
            tickline = replace(&tickline, column, column, &['|']);
            let pad = spaces(column - lines[i].len() as i64);
            lines[i].extend(pad);
            lines[i].extend(chars(&format!("{site} {name}")));
        }
        if lines.is_empty() {
            return;
        }
        let margin = self.margin();
        let tick = format!("{margin}{}\n", text(&tickline));
        let rows: Vec<String> = lines
            .iter()
            .map(|l| format!("{margin}{}\n", text(l)))
            .collect();
        if self.enz_above {
            for row in rows.iter().rev() {
                self.out.insert(row, &[]);
            }
            self.out.insert(&tick, &[]);
        } else {
            self.out.insert(&tick, &[]);
            for row in &rows {
                self.out.insert(row, &[]);
            }
        }
    }

    /// The marks and numbers under or over the DNA: a mark of one kind every so many bases and another kind every so many
    /// more, with the numbers of the second kind written beside them.
    fn index(&mut self, x: usize, linelength: usize) {
        let s = self.s;
        let (sp1, sp2) = (s.index_spacing1 as i64, s.index_spacing2 as i64);
        let start_index = (x + self.start + 1) as i64;
        let first_gap = sp2 - 1 - modulo(start_index - 1, sp2);
        let mut tick = spaces(first_gap);
        tick.push(s.index_char2);
        let unit: Vec<char> = {
            let mut u = spaces(sp2 - 1);
            u.push(s.index_char2);
            u
        };
        for _ in 0..divide(linelength as i64, sp2).max(0) {
            tick.extend(&unit);
        }
        let mut tickline = range(&tick, 0, linelength as i64 - 1);
        let first_index = sp1 - 1 - modulo(start_index - 1, sp1);
        let mut numline = spaces(MARGIN as i64);
        let mut first_num = start_index + sp1 - modulo(start_index, sp1);
        let mut i = first_index;
        while i < linelength as i64 {
            tickline = replace(&tickline, i, i, &[s.index_char1]);
            numline.extend(spaces(sp1));
            let label = chars(&first_num.to_string());
            let keep = MARGIN as i64 - 1 + i - divide(label.len() as i64, 2);
            numline = range(&numline, 0, keep);
            numline.extend(label);
            first_num += sp1;
            i += sp1;
        }
        let tick_row = format!("{}{}\n", self.margin(), text(&tickline));
        let num_row = format!("{}\n", text(&numline));
        if self.index_above {
            if s.index_numbers {
                self.out.insert(&num_row, &[]);
            }
            self.out.insert(&tick_row, &[]);
        } else {
            self.out.insert(&tick_row, &[]);
            if s.index_numbers {
                self.out.insert(&num_row, &[]);
            }
        }
    }

    fn features(&mut self, x: usize, linelength: usize) {
        let first = self.start + x;
        let mut rows = features::Rows::new(&self.tags, self.req.features, self.s, self.row());
        rows.collect(first, linelength, &mut self.cds);
        rows.write(&mut self.out);
    }

    fn finish(mut self) -> TextMap {
        // the tags of the sequence are drawn over the others, in the order the sequence has them
        for tag in &self.tags {
            if self.out.has_tag(&tag.name) {
                self.out.raise(&tag.name);
            }
        }
        let (text, tags) = self.out.finish();
        TextMap {
            text,
            tags,
            dna_rows: self.dna_rows,
        }
    }
}

/// `seq_len_mod` (L58532): a place counted from 1 on a sequence of `len` bases, wrapping round it.
fn seq_len_mod(x: usize, len: usize) -> usize {
    match x % len {
        0 => len,
        rest => rest,
    }
}

/// The base on the other strand: what the second strand line writes under a base (the pairs of `string map` at L19951, which
/// has no entry for u).
fn complement(c: char) -> char {
    match c {
        'A' => 'T',
        'B' => 'V',
        'C' => 'G',
        'D' => 'H',
        'G' => 'C',
        'H' => 'D',
        'K' => 'M',
        'M' => 'K',
        'R' => 'Y',
        'T' => 'A',
        'V' => 'B',
        'Y' => 'R',
        'a' => 't',
        'b' => 'v',
        'c' => 'g',
        'd' => 'h',
        'g' => 'c',
        'h' => 'd',
        'k' => 'm',
        'm' => 'k',
        'r' => 'y',
        't' => 'a',
        'v' => 'b',
        'y' => 'r',
        other => other,
    }
}

/// The translation lines for a whole text (`do_text_map`, L19769-19792). Note that the frames of the six-frame reading with
/// three letters are not moved over to line up, as the one-letter ones are, and that the reverse frames are cut from the end.
fn translations(text: &[u8], code: &[u8], reading: Reading) -> Vec<Vec<char>> {
    let t = |from: usize, forward: bool, three: bool, spacing: usize| -> String {
        genes::translate(
            text.get(from..).unwrap_or(&[]),
            code,
            three,
            forward,
            spacing,
        )
    };
    let one = |from: usize, forward: bool, lead: usize| {
        let mut line = spaces(lead as i64);
        line.extend(chars(&t(from, forward, false, 2)));
        line
    };
    let three = |from: usize| chars(&t(from, true, true, 0));
    let n = text.len();
    // `string range $text 0 end-less`
    let cut = |less: i64| &text[..(n as i64 - less).clamp(0, n as i64) as usize];
    let back = |less: i64| chars(&genes::translate(cut(less), code, true, false, 0));
    match (reading.frames, reading.three_letters) {
        (1, false) => vec![one(0, true, 0)],
        (3, false) => vec![one(0, true, 0), one(1, true, 1), one(2, true, 2)],
        (6, false) => vec![
            one(0, true, 0),
            one(1, true, 1),
            one(2, true, 2),
            one(0, false, 0),
            one(1, false, 1),
            one(2, false, 2),
        ],
        (1, true) => vec![three(0)],
        (3, true) => vec![three(0), spaces_then(1, three(1)), spaces_then(2, three(2))],
        _ => {
            let m = (n % 3) as i64;
            vec![
                three(0),
                three(1),
                three(2),
                back(m),
                back(1 + m),
                back(2 + m),
            ]
        }
    }
}

fn spaces_then(lead: i64, line: Vec<char>) -> Vec<char> {
    let mut out = spaces(lead);
    out.extend(line);
    out
}

/// `linear_enz_sort` (L28745): the sites that fall from `start` to `end` (numbered from 1) with the enzymes' names, by place.
/// A row of one base, where `start` is not before `end`, lists every site of every enzyme: that is how the original goes.
fn linear_enz_sort(
    enzymes: &[EnzymeSites],
    start: usize,
    end: usize,
    seq_length: usize,
) -> Vec<(usize, String)> {
    let end = if end > seq_length {
        end - seq_length
    } else {
        end
    };
    let mut found: Vec<(usize, usize, String)> = Vec::new();
    for e in enzymes {
        if start < end {
            for &site in e.sites {
                if site > end {
                    break;
                }
                if site >= start {
                    found.push((site, site, e.name.to_string()));
                }
            }
        } else {
            for &site in e.sites {
                if site <= end || site >= start {
                    let key = if site <= end { site + seq_length } else { site };
                    found.push((key, site, e.name.to_string()));
                }
            }
        }
    }
    found.sort_by_key(|f| f.0);
    found
        .into_iter()
        .map(|(_, site, name)| (site, name))
        .collect()
}
