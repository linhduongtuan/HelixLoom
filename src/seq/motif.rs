//! Motifs in DNA. A motif is a position weight matrix (PWM): for each of its columns, how likely each base is. It is scanned along
//! sequences on both strands, scored against the background (log-odds, in bits), with the p-value of a score worked out exactly from the
//! distribution of scores of random windows (the way of FIMO, Grant et al. 2011); motifs are found by expectation-maximisation from
//! the words of the sequences, in the simplest form of MEME's method (Bailey and Elkan 1994: at most one site in a sequence); and two
//! motifs are compared by the correlation of their columns. Written from those papers; see `CREDITS.md`.

use super::io::parse_fasta;
use std::collections::{HashMap, HashSet};

pub mod cli;

/// The bases of a column, in the order of its numbers.
pub const BASES: [u8; 4] = *b"ACGT";

/// Scores are worked with in hundredths of a bit when a p-value is found: the error is at most half a hundredth for each column.
const SCALE: f64 = 100.0;
/// How many window-columns the search may go over in all, about five seconds of work: the number of passes over every window of the
/// sequences that it can afford falls as they grow, so that a big set is searched less thoroughly instead of for ever.
const WORK: f64 = 5e8;
/// The most starting words tried, and the least.
const MOST_CANDIDATES: f64 = 3000.0;
const LEAST_CANDIDATES: f64 = 10.0;
/// The fewest iterations of the search from a start, however big the set.
const LEAST_ITERATIONS: usize = 5;
/// Pseudo-counts added to the columns of a motif made from sites, spread by the background.
const PSEUDOCOUNT: f64 = 0.1;

fn code(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Motif {
    pub name: String,
    /// For each column the chance of A, C, G and T, which add up to 1.
    pub probs: Vec<[f64; 4]>,
}

impl Motif {
    pub fn width(&self) -> usize {
        self.probs.len()
    }

    /// The motif of aligned `sites` (letters A, C, G, T, all the same length), with pseudo-counts so that no base is impossible.
    // ponytail: the tests build motifs from sites; a window to type sites in comes with the motif tools that need it
    #[cfg_attr(not(test), allow(dead_code))]
    pub fn from_sites(name: &str, sites: &[&[u8]]) -> Result<Motif, String> {
        let width = sites.first().map_or(0, |s| s.len());
        if width == 0 {
            return Err("There is no site to make a motif of.".into());
        }
        if sites.iter().any(|s| s.len() != width) {
            return Err("The sites are not all the same length.".into());
        }
        let mut counts = vec![[0.0f64; 4]; width];
        for site in sites {
            for (column, &base) in counts.iter_mut().zip(site.iter()) {
                let b = code(base)
                    .ok_or_else(|| format!("\"{}\" is not a base of a site.", base as char))?;
                column[b] += 1.0;
            }
        }
        Ok(Motif {
            name: name.to_string(),
            probs: counts.iter().map(|c| normalise(c, PSEUDOCOUNT)).collect(),
        })
    }

    /// The most likely base of each column.
    pub fn consensus(&self) -> String {
        self.probs
            .iter()
            .map(|c| {
                let best = (0..4).fold(0, |b, k| if c[k] > c[b] { k } else { b });
                BASES[best] as char
            })
            .collect()
    }

    /// How much the motif says, in bits: two for each column that is one base for certain, nothing for a column of equal chances.
    pub fn information(&self) -> f64 {
        self.probs
            .iter()
            .map(|c| {
                2.0 + c
                    .iter()
                    .filter(|&&p| p > 0.0)
                    .map(|&p| p * p.log2())
                    .sum::<f64>()
            })
            .sum()
    }

    /// The motif as MEME's text format reads it (version 4, DNA, both strands), made from `sites` sites: what other programs, such as
    /// FIMO and Tomtom, take.
    pub fn to_meme(&self, background: &[f64; 4], sites: usize) -> String {
        meme_header(background) + &self.meme_block(sites)
    }

    /// The part of a MEME file that is one motif, with the blank line that comes before it.
    fn meme_block(&self, sites: usize) -> String {
        let mut text = format!("\nMOTIF {}\n", self.name.replace(' ', "_"));
        text.push_str(&format!(
            "letter-probability matrix: alength= 4 w= {} nsites= {} E= 0\n",
            self.width(),
            sites
        ));
        for column in &self.probs {
            text.push_str(&format!(
                " {:.6} {:.6} {:.6} {:.6}\n",
                column[0], column[1], column[2], column[3]
            ));
        }
        text
    }

    /// The motif of the opposite strand: columns in the other order, each base replaced by the one that pairs with it.
    pub fn reverse_complement(&self) -> Motif {
        Motif {
            name: self.name.clone(),
            probs: self
                .probs
                .iter()
                .rev()
                .map(|c| [c[3], c[2], c[1], c[0]])
                .collect(),
        }
    }

    /// For each column the log (base 2) of how much likelier each base is in the motif than in the background.
    pub fn log_odds(&self, background: &[f64; 4]) -> Vec<[f64; 4]> {
        self.probs
            .iter()
            .map(|c| std::array::from_fn(|b| (c[b] / background[b]).log2()))
            .collect()
    }
}

/// `counts` as chances, after `pseudo` more counts in all, shared out equally.
fn normalise(counts: &[f64; 4], pseudo: f64) -> [f64; 4] {
    let total: f64 = counts.iter().sum::<f64>() + pseudo;
    std::array::from_fn(|b| (counts[b] + pseudo / 4.0) / total)
}

/// What a MEME file starts with: the version, the alphabet, the strands and the background frequencies.
fn meme_header(background: &[f64; 4]) -> String {
    format!(
        "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA {:.4} C {:.4} G {:.4} T {:.4}\n",
        background[0], background[1], background[2], background[3]
    )
}

/// The motifs found as one MEME file: its header once, then each motif with the number of sites it was made from.
pub fn meme_text(found: &[Found], background: &[f64; 4]) -> String {
    found.iter().fold(meme_header(background), |text, f| {
        text + &f.motif.meme_block(f.sites.len())
    })
}

/// The names and bases of the sequences in `text`: FASTA, or when there is no `>` one sequence for each line. At least two.
pub fn read_sequences(text: &str) -> Result<(Vec<String>, Vec<Vec<u8>>), String> {
    let (names, seqs): (Vec<String>, Vec<Vec<u8>>) = if text.contains('>') {
        parse_fasta(text)
            .into_iter()
            .map(|d| (d.name, d.seq))
            .unzip()
    } else {
        text.lines()
            .map(|l| {
                l.bytes()
                    .filter(u8::is_ascii_alphabetic)
                    .collect::<Vec<u8>>()
            })
            .filter(|s| !s.is_empty())
            .enumerate()
            .map(|(k, s)| (format!("sequence {}", k + 1), s))
            .unzip()
    };
    if seqs.len() < 2 {
        return Err("Give at least two sequences (FASTA, or one on each line).".into());
    }
    Ok((names, seqs))
}

/// The base frequencies of `seqs` (A, C, G, T only; a base never seen still gets a little).
pub fn background(seqs: &[Vec<u8>]) -> [f64; 4] {
    let mut counts = [0.0f64; 4];
    for base in seqs.iter().flatten() {
        if let Some(b) = code(*base) {
            counts[b] += 1.0;
        }
    }
    normalise(&counts, 4.0)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hit {
    /// Where the window starts, counted in the sequence as it is written.
    pub position: usize,
    /// The window matches the motif on the opposite strand.
    pub reverse: bool,
    /// The log-odds score, in bits.
    pub score: f64,
}

/// The windows of `seq` that score at least `min_score` against `motif` (on both strands if asked), in the order of their positions.
/// A window with a letter that is not A, C, G or T is not scored.
pub fn scan(
    seq: &[u8],
    motif: &Motif,
    background: &[f64; 4],
    min_score: f64,
    both_strands: bool,
) -> Vec<Hit> {
    let width = motif.width();
    let forward = motif.log_odds(background);
    let reverse = motif.reverse_complement().log_odds(background);
    if width == 0 || seq.len() < width {
        return Vec::new();
    }
    let codes: Vec<Option<usize>> = seq.iter().map(|&b| code(b)).collect();
    let mut hits = Vec::new();
    for position in 0..=seq.len() - width {
        let window = &codes[position..position + width];
        if window.iter().any(Option::is_none) {
            continue;
        }
        let score = |matrix: &[[f64; 4]]| -> f64 {
            window
                .iter()
                .zip(matrix)
                .map(|(b, column)| column[b.unwrap_or(0)])
                .sum()
        };
        let f = score(&forward);
        if f >= min_score {
            hits.push(Hit {
                position,
                reverse: false,
                score: f,
            });
        }
        if both_strands {
            let r = score(&reverse);
            if r >= min_score {
                hits.push(Hit {
                    position,
                    reverse: true,
                    score: r,
                });
            }
        }
    }
    hits
}

/// The distribution of the score of a window of random bases (drawn with the background's frequencies) against a motif, with the
/// scores rounded to hundredths of a bit, from which the p-value of a score is read exactly.
pub struct ScoreDistribution {
    /// The scaled score of the lowest bin.
    lowest: i64,
    /// For each bin, the chance of a score at least that high.
    at_least: Vec<f64>,
}

impl ScoreDistribution {
    pub fn new(motif: &Motif, background: &[f64; 4]) -> ScoreDistribution {
        let columns: Vec<[i64; 4]> = motif
            .log_odds(background)
            .iter()
            .map(|c| std::array::from_fn(|b| (c[b] * SCALE).round() as i64))
            .collect();
        let lowest: i64 = columns.iter().map(|c| *c.iter().min().unwrap_or(&0)).sum();
        let highest: i64 = columns.iter().map(|c| *c.iter().max().unwrap_or(&0)).sum();
        let size = (highest - lowest) as usize + 1;
        // after each column: the chance of every total so far, counted from the lowest possible
        let mut chance = vec![0.0f64; size];
        chance[0] = 1.0;
        let mut reached = 0usize;
        let mut floor = 0i64;
        for column in &columns {
            let least = *column.iter().min().unwrap_or(&0);
            let most = *column.iter().max().unwrap_or(&0);
            let mut next = vec![0.0f64; size];
            for (at, &p) in chance.iter().enumerate().take(reached + 1) {
                if p == 0.0 {
                    continue;
                }
                for (b, &s) in column.iter().enumerate() {
                    // the new total, counted from the new lowest (the floor moves up by this column's least)
                    let to = at as i64 + (s - least);
                    next[to as usize] += p * background[b];
                }
            }
            reached += (most - least) as usize;
            floor += least;
            chance = next;
        }
        debug_assert_eq!(floor, lowest);
        let mut at_least = chance;
        for k in (0..at_least.len().saturating_sub(1)).rev() {
            at_least[k] += at_least[k + 1];
        }
        ScoreDistribution { lowest, at_least }
    }

    /// The chance that a random window scores at least `score` (to a hundredth of a bit for each column).
    pub fn p_value(&self, score: f64) -> f64 {
        // a score that is exactly some hundredths must not be rounded up past them by the error of the division that made it
        let scaled = (score * SCALE - 1e-6).ceil() as i64 - self.lowest;
        if scaled <= 0 {
            1.0
        } else if scaled as usize >= self.at_least.len() {
            0.0
        } else {
            self.at_least[scaled as usize].min(1.0)
        }
    }

    /// The lowest score whose p-value is at most `p`.
    pub fn score_for(&self, p: f64) -> f64 {
        let at = self
            .at_least
            .iter()
            .position(|&chance| chance <= p)
            .unwrap_or(self.at_least.len());
        (self.lowest + at as i64) as f64 / SCALE
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Settings {
    pub width: usize,
    /// How many motifs to look for, one after the other (the places of a motif found are hidden from the next search).
    pub motifs: usize,
    pub both_strands: bool,
    /// How many starting words are taken through the whole of the search.
    pub seeds: usize,
    pub iterations: usize,
    /// How much work the search may spend, as a multiple of the usual budget of about five seconds on a big set (a small set is
    /// searched in full whatever this is).
    pub effort: f64,
}

/// The least and the most effort a search may be given.
pub const EFFORT_RANGE: (f64, f64) = (0.01, 1000.0);

/// The most memory the windows of the sequences may take: every window of every strand is held as the numbers of its bases.
const MAX_WINDOW_BYTES: u128 = 4 << 30;

/// Refuse a search whose windows would not fit in memory (`bases` bases, motifs `width` columns wide, on both strands or one),
/// before any is made. A window costs about 64 bytes and 8 for each column.
pub fn check_size(bases: usize, width: usize, both_strands: bool) -> Result<(), String> {
    let windows = bases as u128 * if both_strands { 2 } else { 1 };
    let bytes = windows * (64 + 8 * width as u128);
    if bytes > MAX_WINDOW_BYTES {
        return Err(format!(
            "{bases} bases need about {} GB to hold the windows of a motif {width} columns wide; at most {} GB are taken: use fewer bases or a narrower motif.",
            bytes >> 30,
            MAX_WINDOW_BYTES >> 30
        ));
    }
    Ok(())
}

impl Default for Settings {
    fn default() -> Settings {
        Settings {
            width: 8,
            motifs: 1,
            both_strands: true,
            seeds: 12,
            iterations: 60,
            effort: 1.0,
        }
    }
}

/// One place a motif was found.
#[derive(Clone, Debug, PartialEq)]
pub struct Site {
    pub sequence: usize,
    pub position: usize,
    pub reverse: bool,
    /// The log-odds score of the window, in bits.
    pub score: f64,
    /// How likely it is, after the search, that this is a site of the motif.
    pub probability: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Found {
    pub motif: Motif,
    /// In each sequence the best window, where it is more likely than not a site.
    pub sites: Vec<Site>,
    /// How much likelier the sequences are with the motif than with the background alone, in bits: the larger, the better the motif.
    pub gain: f64,
}

/// Every window of the sequence as the numbers of its bases, with the strand it is read from and where it starts; none for a window
/// with a letter that is not a base.
fn windows(
    seq: &[Vec<Option<usize>>],
    width: usize,
    both: bool,
) -> Vec<Vec<(usize, bool, Vec<usize>)>> {
    seq.iter()
        .map(|codes| {
            let mut out = Vec::new();
            if codes.len() < width {
                return out;
            }
            for position in 0..=codes.len() - width {
                let Some(word) = codes[position..position + width]
                    .iter()
                    .copied()
                    .collect::<Option<Vec<usize>>>()
                else {
                    continue;
                };
                if both {
                    out.push((position, true, word.iter().rev().map(|&b| 3 - b).collect()));
                }
                out.push((position, false, word));
            }
            out
        })
        .collect()
}

/// One round of expectation-maximisation: the chance of a site at each window (at most one site in a sequence), and from those the
/// new motif and the new chance that a sequence has a site. Returns the gain of the motif it started from, in bits.
fn step(
    words: &[Vec<(usize, bool, Vec<usize>)>],
    probs: &[[f64; 4]],
    background: &[f64; 4],
    prior: f64,
) -> (Vec<[f64; 4]>, f64, f64, Vec<Vec<f64>>) {
    let mut counts = vec![[0.0f64; 4]; probs.len()];
    let (mut total_site, mut sequences, mut gain) = (0.0f64, 0.0f64, 0.0f64);
    let mut z_all = Vec::with_capacity(words.len());
    for sequence in words {
        if sequence.is_empty() {
            z_all.push(Vec::new());
            continue;
        }
        sequences += 1.0;
        let ratio: Vec<f64> = sequence
            .iter()
            .map(|(_, _, word)| {
                word.iter()
                    .zip(probs)
                    .map(|(&b, column)| (column[b] / background[b]).ln())
                    .sum::<f64>()
                    .exp()
            })
            .collect();
        // the prior chance of a site is shared among the windows of the sequence; both strands of a place are two windows
        let each = prior / sequence.len() as f64;
        let evidence = (1.0 - prior) + each * ratio.iter().sum::<f64>();
        gain += evidence.log2();
        let z: Vec<f64> = ratio.iter().map(|r| each * r / evidence).collect();
        for ((_, _, word), &weight) in sequence.iter().zip(&z) {
            for (k, &b) in word.iter().enumerate() {
                counts[k][b] += weight;
            }
        }
        total_site += z.iter().sum::<f64>();
        z_all.push(z);
    }
    let next: Vec<[f64; 4]> = counts
        .iter()
        .map(|c| normalise(c, PSEUDOCOUNT * 4.0))
        .collect();
    (
        next,
        (total_site / sequences.max(1.0)).clamp(0.01, 0.99),
        gain,
        z_all,
    )
}

/// What the search can afford over `cost` window-columns (windows times width) for each pass over all of them: how many starting words to
/// try, how many starts to follow, and how many iterations each. The settings are what a small set gets whole.
pub(crate) fn plan(cost: f64, settings: &Settings) -> (usize, usize, usize) {
    let passes = (WORK * settings.effort / cost.max(1.0)).max(1.0);
    let candidates = (passes * 0.4).clamp(LEAST_CANDIDATES, MOST_CANDIDATES) as usize;
    let left = (passes * 0.6).max(LEAST_ITERATIONS as f64);
    let iterations = settings.iterations.min(left as usize).max(LEAST_ITERATIONS);
    let seeds = settings
        .seeds
        .max(1)
        .min(((left / iterations as f64).floor() as usize).max(1));
    (candidates, seeds, iterations)
}

/// A number for a word that is the same in every run (the order of equal candidates must not depend on a hash seed).
fn word_hash(word: &[usize]) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    word.hash(&mut hasher);
    hasher.finish()
}

/// A starting motif from a word: the word's base is likely, the others share the rest.
fn seed_motif(word: &[usize]) -> Vec<[f64; 4]> {
    word.iter()
        .map(|&b| std::array::from_fn(|k| if k == b { 0.5 } else { 0.5 / 3.0 }))
        .collect()
}

/// Look for `settings.motifs` motifs of `settings.width` columns in `seqs`. `Err` when no sequence is long enough or has a base.
pub fn discover(seqs: &[Vec<u8>], settings: &Settings) -> Result<Vec<Found>, String> {
    let width = settings.width;
    if !(3..=40).contains(&width) {
        return Err("A motif is between 3 and 40 columns wide.".into());
    }
    if !(EFFORT_RANGE.0..=EFFORT_RANGE.1).contains(&settings.effort) {
        return Err(format!(
            "The effort of the search is from {} to {}, not {}.",
            EFFORT_RANGE.0, EFFORT_RANGE.1, settings.effort
        ));
    }
    check_size(
        seqs.iter().map(Vec::len).sum(),
        width,
        settings.both_strands,
    )?;
    if seqs.iter().all(|s| s.len() < width) {
        return Err(format!(
            "No sequence is as long as the motif ({width} bases)."
        ));
    }
    let bg = background(seqs);
    let mut working: Vec<Vec<Option<usize>>> = seqs
        .iter()
        .map(|s| s.iter().map(|&b| code(b)).collect())
        .collect();
    let mut found = Vec::new();
    for round in 0..settings.motifs.max(1) {
        let words = windows(&working, width, settings.both_strands);
        // the starting words: the distinct words of the sequences, each judged by how well one step of the search from it does
        // The starting words: every different word of the sequences, with the number of sequences that hold it. A word that many hold is
        // likelier to be (a copy of) a site than one that one holds, so the words held by the most come first, and only as many are tried
        // as the work allows (ties go in an order that is the same every time but is not the order of the sequences).
        let mut holders: HashMap<&Vec<usize>, usize> = HashMap::new();
        for sequence in &words {
            let held: HashSet<&Vec<usize>> = sequence.iter().map(|(_, _, word)| word).collect();
            for word in held {
                *holders.entry(word).or_insert(0) += 1;
            }
        }
        if holders.is_empty() {
            break;
        }
        let cost = words.iter().map(Vec::len).sum::<usize>() as f64 * width as f64;
        let (most_candidates, seeds, iterations) = plan(cost, settings);
        let mut candidates: Vec<(usize, u64, &Vec<usize>)> = holders
            .into_iter()
            .map(|(word, count)| (count, word_hash(word), word))
            .collect();
        candidates.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
        let mut scored: Vec<(f64, &Vec<usize>)> = candidates
            .iter()
            .take(most_candidates)
            .map(|(_, _, word)| (step(&words, &seed_motif(word), &bg, 0.5).2, *word))
            .collect();
        scored.sort_by(|a, b| b.0.total_cmp(&a.0));
        let mut best: Option<(f64, Vec<[f64; 4]>, f64)> = None;
        for (_, word) in scored.iter().take(seeds) {
            let (mut probs, mut prior, mut gain) = (seed_motif(word), 0.5, f64::MIN);
            for _ in 0..iterations {
                let (next, next_prior, g, _) = step(&words, &probs, &bg, prior);
                let change = next
                    .iter()
                    .zip(&probs)
                    .map(|(a, b)| (0..4).map(|k| (a[k] - b[k]).abs()).fold(0.0, f64::max))
                    .fold(0.0, f64::max);
                (probs, prior, gain) = (next, next_prior, g);
                if change < 1e-5 {
                    break;
                }
            }
            if best.as_ref().is_none_or(|b| gain > b.0) {
                best = Some((gain, probs, prior));
            }
        }
        let Some((_, probs, prior)) = best else { break };
        let (_, _, gain, z) = step(&words, &probs, &bg, prior);
        let motif = Motif {
            name: format!("motif {}", round + 1),
            probs,
        };
        let mut sites = Vec::new();
        for (s, (sequence, zs)) in words.iter().zip(&z).enumerate() {
            let Some((k, &p)) = zs.iter().enumerate().max_by(|a, b| a.1.total_cmp(b.1)) else {
                continue;
            };
            if p < 0.5 {
                continue;
            }
            let (position, reverse, word) = &sequence[k];
            let score: f64 = word
                .iter()
                .zip(&motif.probs)
                .map(|(&b, c)| (c[b] / bg[b]).log2())
                .sum();
            sites.push(Site {
                sequence: s,
                position: *position,
                reverse: *reverse,
                score,
                probability: p,
            });
            // hide the place from the next search
            for cell in working[s][*position..*position + width].iter_mut() {
                *cell = None;
            }
        }
        found.push(Found { motif, sites, gain });
    }
    Ok(found)
}

/// How alike two motifs are: the best mean Pearson correlation of their columns over the ways they can overlap by at least
/// `least` columns (or all of the narrower), on either strand of the second. 1 for the same motif, near 0 for unrelated ones.
pub fn compare(a: &Motif, b: &Motif, least: usize) -> f64 {
    let column_correlation = |x: &[f64; 4], y: &[f64; 4]| -> f64 {
        let (mx, my) = (0.25, 0.25);
        let (mut sxy, mut sxx, mut syy) = (0.0, 0.0, 0.0);
        for k in 0..4 {
            sxy += (x[k] - mx) * (y[k] - my);
            sxx += (x[k] - mx).powi(2);
            syy += (y[k] - my).powi(2);
        }
        if sxx == 0.0 || syy == 0.0 {
            0.0
        } else {
            sxy / (sxx * syy).sqrt()
        }
    };
    let least = least.max(1).min(a.width()).min(b.width());
    let mut best = -1.0f64;
    for other in [b.clone(), b.reverse_complement()] {
        let (wa, wb) = (a.width() as i64, other.width() as i64);
        for shift in -(wb - least as i64)..=(wa - least as i64) {
            // column k of `a` against column k - shift of `other`
            let pairs: Vec<(usize, usize)> = (0..wa)
                .filter(|k| (0..wb).contains(&(k - shift)))
                .map(|k| (k as usize, (k - shift) as usize))
                .collect();
            if pairs.len() < least {
                continue;
            }
            let mean = pairs
                .iter()
                .map(|&(i, j)| column_correlation(&a.probs[i], &other.probs[j]))
                .sum::<f64>()
                / pairs.len() as f64;
            best = best.max(mean);
        }
    }
    best
}

#[cfg(test)]
mod tests;
