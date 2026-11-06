//! Built-in demo data: a synthetic circular plasmid and sequencing-style reads sampled from it.
//! The sequence is generated, not taken from any real plasmid.

use crate::seq::assemble::Read;
use crate::seq::dna::{find_all, revcom, FindOpts};
use crate::seq::doc::Doc;
use crate::seq::features::Feature;

/// Multiple cloning site: EcoRI SacI KpnI SmaI BamHI XbaI SalI PstI SphI HindIII.
const MCS: &[u8] = b"GAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTT";
const MCS_ENZYME_SITES: &[&[u8]] = &[
    b"GAATTC", b"GAGCTC", b"GGTACC", b"CCCGGG", b"GGATCC", b"TCTAGA", b"GTCGAC", b"CTGCAG",
    b"GCATGC", b"AAGCTT",
];

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }

    fn dna(&mut self, len: usize) -> Vec<u8> {
        (0..len)
            .map(|_| b"ACGT"[(self.next() % 4) as usize])
            .collect()
    }

    /// `codons` random sense codons (no stop) between ATG and a stop codon.
    fn orf(&mut self, codons: usize) -> Vec<u8> {
        let mut out = b"ATG".to_vec();
        while out.len() < 3 + codons * 3 {
            let c = self.dna(3);
            if !matches!(c.as_slice(), b"TAA" | b"TAG" | b"TGA") {
                out.extend(c);
            }
        }
        out.extend_from_slice(b"TAA");
        out
    }
}

/// Make sure the MCS enzymes cut only inside the MCS by nudging accidental copies.
fn scrub(seq: &mut [u8], protected: (usize, usize)) {
    let opts = FindOpts {
        both_strands: true,
        ..Default::default()
    };
    for site in MCS_ENZYME_SITES {
        while let Some(hit) = find_all(seq, site, &opts)
            .into_iter()
            .find(|h| h.start < protected.0 || h.start >= protected.1)
        {
            let mid = hit.start + hit.len / 2;
            seq[mid] = if seq[mid] == b'A' { b'C' } else { b'A' };
        }
    }
}

pub fn demo_plasmid() -> Doc {
    let mut rng = Lcg(2026);
    let mut seq: Vec<u8> = Vec::new();
    let mut features: Vec<Feature> = Vec::new();
    let mut add = |seq: &mut Vec<u8>, part: &[u8], name: &str, kind: &str, reverse: bool| {
        let start = seq.len();
        seq.extend_from_slice(part);
        if !name.is_empty() {
            features.push(Feature::new(name, kind, start, seq.len(), reverse));
        }
    };
    add(&mut seq, &rng.dna(150), "", "", false);
    add(
        &mut seq,
        b"TAATACGACTCACTATAGGG",
        "T7 promoter",
        "promoter",
        false,
    );
    add(&mut seq, &rng.dna(14), "", "", false);
    add(
        &mut seq,
        b"AATTGTGAGCGGATAACAATT",
        "lac operator",
        "misc_binding",
        false,
    );
    add(&mut seq, &rng.dna(18), "", "", false);
    let mcs_start = seq.len();
    add(&mut seq, MCS, "MCS", "misc_feature", false);
    let mcs_end = seq.len();
    add(&mut seq, &rng.dna(12), "", "", false);
    add(
        &mut seq,
        b"GACTACAAGGACGACGATGACAAG",
        "FLAG",
        "misc_feature",
        false,
    );
    add(&mut seq, &rng.dna(10), "", "", false);
    add(&mut seq, &rng.orf(190), "Reporter ORF (demo)", "CDS", false);
    add(&mut seq, &rng.dna(30), "", "", false);
    add(
        &mut seq,
        b"CACCACCACCACCACCAC",
        "6xHis",
        "misc_feature",
        false,
    );
    add(&mut seq, &rng.dna(26), "", "", false);
    add(
        &mut seq,
        b"CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG",
        "T7 terminator",
        "terminator",
        false,
    );
    add(&mut seq, &rng.dna(120), "", "", false);
    add(
        &mut seq,
        b"ATAACTTCGTATAATGTATGCTATACGAAGTTAT",
        "loxP",
        "misc_recomb",
        false,
    );
    add(&mut seq, &rng.dna(90), "", "", false);
    add(
        &mut seq,
        &revcom(&rng.orf(250)),
        "Marker ORF (demo)",
        "CDS",
        true,
    );
    add(&mut seq, &rng.dna(80), "", "", false);
    add(
        &mut seq,
        &rng.dna(560),
        "Origin (demo)",
        "rep_origin",
        false,
    );
    add(&mut seq, &rng.dna(60), "", "", false);
    scrub(&mut seq, (mcs_start, mcs_end));
    let mut doc = Doc::new("pDemo", seq);
    doc.circular = true;
    doc.features = features;
    doc.definition = "Synthetic demo plasmid generated for HelixLoom.".into();
    doc.dirty = false;
    doc
}

/// 300 synthetic Illumina-style reads from the demo plasmid for the read-cleaning window: 100
/// bases with quality falling towards the 3' end, every seventh read a short insert followed by
/// the TruSeq adapter, and a few with runs of N.
pub fn demo_fastq() -> String {
    const ADAPTER: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    let plasmid = demo_plasmid().seq;
    let mut state = 12345u32;
    let mut next = move |n: u32| {
        state = state.wrapping_mul(1_103_515_245).wrapping_add(12345);
        (state >> 16) % n
    };
    (0..300)
        .map(|i| {
            let start = next(plasmid.len() as u32 - 150) as usize;
            let insert = if i % 7 == 0 {
                40 + next(40) as usize
            } else {
                100
            };
            let mut seq = plasmid[start..start + insert].to_vec();
            seq.extend(ADAPTER.iter().take(100 - insert));
            // past the adapter a read carries whatever follows it
            while seq.len() < 100 {
                seq.push(b"ACGT"[next(4) as usize]);
            }
            let mut qual: Vec<u8> = (0..100u32)
                .map(|p| (38 + next(5) as i32 - 2 - (p * 26 / 100) as i32).clamp(2, 41) as u8)
                .collect();
            if i % 23 == 0 {
                for p in 88..95 {
                    (seq[p], qual[p]) = (b'N', 2);
                }
            }
            let text: String = qual.iter().map(|&q| char::from(q + 33)).collect();
            format!("@demo_{i}\n{}\n+\n{text}\n", String::from_utf8_lossy(&seq))
        })
        .collect()
}

/// Six related sequences for the tree view: the demo plasmid's first 600 bases with growing
/// numbers of substitutions (the second and third share their first changes; the fourth and sixth
/// branch off on their own).
pub fn tree_demo() -> Vec<Doc> {
    let base: Vec<u8> = demo_plasmid().seq[..600].to_vec();
    let variants = [
        ("Isolate 1", 0, 0),
        ("Isolate 2", 4, 0),
        ("Isolate 3", 8, 0),
        ("Isolate 4", 8, 100),
        ("Isolate 5", 20, 0),
        ("Isolate 6", 40, 200),
    ];
    variants
        .into_iter()
        .map(|(name, changes, offset)| {
            let mut seq = base.clone();
            for i in 0..changes {
                // 41 shares no factor with 600, so the changed positions never repeat
                let at = (i * 41 + offset + 7) % seq.len();
                seq[at] = match seq[at] {
                    b'A' => b'G',
                    b'G' => b'A',
                    b'C' => b'T',
                    _ => b'C',
                };
            }
            Doc::new(name, seq)
        })
        .collect()
}

/// A made-up tree of 24 strains in three groups (East_, West_, North_), with branch lengths and support
/// values, for the tree viewer's example.
pub fn demo_tree() -> crate::seq::tree::Node {
    use crate::seq::tree::Node;
    fn unit(rng: &mut Lcg) -> f64 {
        (rng.next() % 1000) as f64 / 1000.0
    }
    /// A random branching of `n` tips named `prefix_01`, `prefix_02`, ..., deeper branches shorter.
    fn clade(rng: &mut Lcg, prefix: &str, n: usize, counter: &mut usize, scale: f64) -> Node {
        let length = scale * (0.02 + 0.04 * unit(rng));
        if n == 1 {
            *counter += 1;
            return Node {
                name: format!("{prefix}_{:02}", *counter),
                length,
                support: None,
                children: Vec::new(),
            };
        }
        let left = 1 + (rng.next() as usize) % (n - 1);
        let children = vec![
            clade(rng, prefix, left, counter, scale * 0.8),
            clade(rng, prefix, n - left, counter, scale * 0.8),
        ];
        Node {
            name: String::new(),
            length,
            support: Some(60.0 + (rng.next() % 41) as f64),
            children,
        }
    }
    let mut rng = Lcg(2024);
    let join = |children: Vec<Node>, length: f64, support: f64| Node {
        name: String::new(),
        length,
        support: Some(support),
        children,
    };
    let east = clade(&mut rng, "East", 9, &mut 0, 1.0);
    let west = clade(&mut rng, "West", 8, &mut 0, 1.0);
    let north = clade(&mut rng, "North", 7, &mut 0, 1.0);
    let inner = join(vec![east, west], 0.05, 100.0);
    let mut root = join(vec![inner, north], 0.0, 100.0);
    root.support = None;
    root
}

/// A made-up Sanger trace of the demo plasmid's bases 151 to 630 for the Signals page's example: peaks a dozen
/// scans apart, a messy start, a clean middle and a weakening end, with quality values to match. The signal is
/// computed, not recorded.
pub fn demo_trace() -> crate::seq::abif::Trace {
    const FIRST: usize = 150;
    const LAST: usize = 630;
    const ORDER: [u8; 4] = *b"GATC";
    const WIDTH: f32 = 2.6;
    let bases: Vec<u8> = demo_plasmid().seq[FIRST..LAST].to_vec();
    let n = bases.len();
    let mut rng = Lcg(5);
    let mut at = 14.0f32;
    let peaks: Vec<u32> = (0..n)
        .map(|_| {
            let here = at.round() as u32;
            at += 11.0 + (rng.next() % 40) as f32 / 10.0;
            here
        })
        .collect();
    let scans = peaks.last().map_or(0, |&p| p as usize) + 24;
    let mut signal: [Vec<f32>; 4] = std::array::from_fn(|_| vec![0.0; scans]);
    let mut quality = Vec::with_capacity(n);
    for (i, (&base, &peak)) in bases.iter().zip(&peaks).enumerate() {
        // weak and ragged for the first bases, strong in the middle, fading towards the end
        let health = if i < 60 {
            0.05 + 0.95 * i as f32 / 60.0
        } else if i + 110 > n {
            (1.0 - (i + 110 - n) as f32 / 110.0).max(0.1)
        } else {
            1.0
        };
        let height = health * (650.0 + (rng.next() % 450) as f32);
        let k = ORDER.iter().position(|&o| o == base).unwrap_or(0);
        let peak = peak as usize;
        let from = peak.saturating_sub(9);
        for (s, level) in signal[k][from..(peak + 10).min(scans)]
            .iter_mut()
            .enumerate()
        {
            let d = ((from + s) as f32 - peak as f32) / WIDTH;
            *level += height * (-d * d / 2.0).exp();
        }
        // the neighbouring base's channel catches a little of every peak, more when the run is poor
        signal[(k + 1) % 4][peak] += height * 0.12 * (1.5 - health);
        quality.push((health * 52.0 + (rng.next() % 9) as f32).clamp(3.0, 60.0) as u8);
    }
    let channels = signal.map(|c| {
        c.iter()
            .map(|&v| (v + (rng.next() % 13) as f32 - 6.0) as i16)
            .collect()
    });
    crate::seq::abif::Trace {
        sample: "pDemo_T7_forward (example)".into(),
        instrument: "synthetic".into(),
        order: ORDER,
        channels,
        bases,
        quality,
        peaks,
        ..crate::seq::abif::Trace::default()
    }
}

/// Reads sampled around the demo plasmid: both strands, a few substitutions, one short contained
/// read and one unrelated read.
pub fn demo_reads() -> Vec<Read> {
    let plasmid = demo_plasmid().seq;
    let n = plasmid.len();
    let mut rng = Lcg(77);
    let mut reads = Vec::new();
    let (width, step) = (760usize, 230usize);
    let mut start = 40usize;
    for i in 0..n.div_ceil(step) + 1 {
        let w = width - (i * 37) % 120;
        let mut s: Vec<u8> = (0..w).map(|k| plasmid[(start + k) % n]).collect();
        if i % 4 == 2 {
            let at = 100 + (i * 53) % 300;
            s[at] = if s[at] == b'G' { b'T' } else { b'G' };
        }
        let reverse = i % 2 == 1;
        let name = format!("read{:02}_{}", i + 1, if reverse { "R" } else { "F" });
        reads.push(Read {
            name,
            seq: if reverse { revcom(&s) } else { s },
        });
        start += step;
    }
    let inner: Vec<u8> = (0..180).map(|k| plasmid[(900 + k) % n]).collect();
    reads.push(Read {
        name: "short_internal_F".into(),
        seq: inner,
    });
    reads.push(Read {
        name: "unrelated_contaminant".into(),
        seq: rng.dna(420),
    });
    reads
}

/// A made-up binary tree of `tips` tips (`t001`, `t002`, …) with branch lengths, as Newick text: the PhyloCNN window's
/// example where the Zurich HIV tree is not installed (that tree belongs to the authors of PhyloCNN and is not part of
/// this program). Deeper branches are shorter, as in a tree of an outbreak.
pub fn made_up_tree(tips: usize) -> String {
    /// A random branching of `n` tips; the part of the text after the brackets is the branch to it (none for the root).
    fn clade(rng: &mut Lcg, n: usize, counter: &mut usize, scale: f64, root: bool) -> String {
        let length = if root {
            String::new()
        } else {
            format!(
                ":{:.4}",
                scale * (0.2 + 0.8 * (rng.next() % 1000) as f64 / 1000.0)
            )
        };
        if n <= 1 {
            *counter += 1;
            return format!("t{:03}{length}", *counter);
        }
        let left = 1 + (rng.next() as usize) % (n - 1);
        let first = clade(rng, left, counter, scale * 0.97, false);
        let second = clade(rng, n - left, counter, scale * 0.97, false);
        format!("({first},{second}){length}")
    }
    let mut rng = Lcg(1_987);
    format!("{};", clade(&mut rng, tips.max(2), &mut 0, 1.0, true))
}

/// `count` copies of the example tree of the PhyloCNN window with their branch lengths shaken a little, as Newick text
/// one tree per line: a made-up posterior sample for the window's summary of many trees.
pub fn tree_sample(count: usize) -> String {
    let (example, _) = crate::seq::phylocnn::example_tree();
    let Ok(mut parsed) = crate::seq::newick::NwTree::parse_all(&example) else {
        return String::new();
    };
    let tree = parsed.remove(0);
    (0..count)
        .map(|k| {
            let mut copy = tree.clone();
            for (i, node) in copy.nodes.iter_mut().enumerate() {
                node.dist *= 1.0 + 0.1 * ((i + 7 * k) as f64).sin();
            }
            copy.newick()
        })
        .collect::<Vec<_>>()
        .join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::assemble::{assemble, Params};
    use crate::seq::enzymes::{all_sites, default_enzymes};

    #[test]
    fn the_demo_plasmid_is_circular_with_features_and_unique_mcs_sites() {
        let d = demo_plasmid();
        assert!(d.circular);
        assert!((2000..3000).contains(&d.len()), "length {}", d.len());
        assert!(d.features.len() >= 10);
        let enzymes = default_enzymes();
        let sites = all_sites(&enzymes, &d.seq, true, false);
        for name in [
            "EcoRI", "SacI", "KpnI", "SmaI", "BamHI", "XbaI", "SalI", "PstI", "SphI", "HindIII",
        ] {
            let i = enzymes.iter().position(|e| e.name == name).unwrap();
            assert_eq!(sites[i].len(), 1, "{name} should cut once");
        }
        for f in &d.features {
            assert!(f.end() <= d.len(), "{} runs past the end", f.name);
        }
    }

    #[test]
    fn the_demo_orfs_are_real_open_reading_frames() {
        let d = demo_plasmid();
        for f in d.features.iter().filter(|f| f.kind == "CDS") {
            let bases = d.feature_bases(f);
            let protein = crate::seq::dna::translate(&bases, crate::seq::dna::STANDARD_CODE);
            assert_eq!(*protein.last().unwrap(), b'*');
            assert_eq!(
                protein.iter().filter(|&&a| a == b'*').count(),
                1,
                "{}",
                f.name
            );
        }
    }

    #[test]
    fn the_made_up_tree_has_the_tips_asked_for_and_is_one_the_networks_can_encode() {
        let text = made_up_tree(200);
        let trees = crate::seq::newick::NwTree::parse_all(&text).unwrap();
        assert_eq!(trees.len(), 1);
        assert_eq!(trees[0].leaf_count(), 200);
        // two tips per inner node less one: 200 tips and 199 inner nodes
        assert_eq!(
            crate::seq::phylocnn::encode(&trees[0]).unwrap().rows.len(),
            399
        );
        assert_eq!(text, made_up_tree(200), "the same tree every time");
        assert_eq!(
            crate::seq::newick::NwTree::parse_all(&made_up_tree(1))
                .unwrap()
                .remove(0)
                .leaf_count(),
            2,
            "a tree has at least two tips"
        );
    }

    #[test]
    fn the_tree_sample_has_one_usable_tree_per_line() {
        let text = tree_sample(5);
        let trees = crate::seq::newick::NwTree::parse_all(&text).unwrap();
        assert_eq!(trees.len(), 5);
        assert!(trees.iter().all(|t| t.leaf_count() == 200));
        assert_ne!(text.lines().next(), text.lines().nth(1));
    }

    #[test]
    fn the_demo_tree_has_three_groups_of_named_tips_and_survives_newick() {
        let tree = demo_tree();
        assert_eq!(tree.leaf_count(), 24);
        let parsed = crate::seq::newick::NwTree::parse_all(&tree.newick())
            .unwrap()
            .remove(0);
        assert_eq!(parsed.leaf_count(), 24);
        let names: Vec<String> = parsed
            .nodes
            .iter()
            .filter(|n| n.children.is_empty())
            .map(|n| n.name.clone())
            .collect();
        for group in ["East_", "West_", "North_"] {
            assert!(names.iter().any(|n| n.starts_with(group)), "{group}");
        }
        // the same tree every time
        assert_eq!(demo_tree().newick(), tree.newick());
    }

    #[test]
    fn the_demo_trace_has_a_peak_in_the_right_colour_under_every_base_and_weak_ends() {
        let t = demo_trace();
        assert_eq!(t.bases.len(), 480);
        assert_eq!((t.quality.len(), t.peaks.len()), (480, 480));
        assert!(
            t.peaks.windows(2).all(|w| w[0] < w[1]),
            "peaks run in order"
        );
        let right = t
            .bases
            .iter()
            .zip(&t.peaks)
            .filter(|&(&base, &peak)| {
                let strongest = (0..4)
                    .max_by_key(|&k| t.channels[k][peak as usize])
                    .expect("four channels");
                t.channel_of(base) == Some(strongest)
            })
            .count();
        assert!(
            right >= 475,
            "{right} of 480 peaks carry their base's colour"
        );
        // weak start and end, strong middle: trimming at the default cutoff keeps the middle only
        let (from, to) = t.good_stretch(0.05).expect("a good stretch");
        assert!(
            (4..40).contains(&from) && (380..478).contains(&to),
            "{from}..{to}"
        );
        assert!(t.quality[240] > 40);
        // the example is the plasmid's own bases
        assert_eq!(t.bases, demo_plasmid().seq[150..630].to_vec());
    }

    #[test]
    fn the_demo_reads_assemble_back_into_the_plasmid() {
        let reads = demo_reads();
        let asm = assemble(&reads, &Params::reads());
        assert_eq!(asm.contigs.len(), 1, "contigs: {}", asm.contigs.len());
        let c = &asm.contigs[0];
        assert!(c.circular);
        let plasmid = demo_plasmid().seq;
        assert_eq!(c.consensus.len(), plasmid.len());
        let doubled: Vec<u8> = plasmid.iter().chain(plasmid.iter()).copied().collect();
        let rc_doubled = revcom(&doubled);
        let matches = |hay: &[u8]| {
            hay.windows(plasmid.len())
                .any(|w| w == c.consensus.as_slice())
        };
        assert!(matches(&doubled) || matches(&rc_doubled));
        assert_eq!(asm.unassembled, vec![reads.len() - 1]);
    }
}
