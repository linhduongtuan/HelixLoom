//! The workflow's steps. Each uses the installed program when there is one and the program's own built-in
//! engine otherwise: MAFFT or the centre-star alignment, IQ-TREE or the likelihood search, primer3 or the
//! port of it.

use super::tools::{self, Installed};
use crate::seq::msa::{center_star, Msa};
use crate::seq::newick::NwTree;
use crate::seq::primer3::Settings;
use std::path::Path;
use web_time::Instant;

pub type Named = (String, Vec<u8>);

/// Names that are safe in every program's files: s1, s2, ... The real names are put back afterwards.
fn ids(n: usize) -> Vec<String> {
    (1..=n).map(|i| format!("s{i}")).collect()
}

pub fn fasta_text(seqs: &[(String, Vec<u8>)]) -> String {
    seqs.iter()
        .map(|(name, bases)| format!(">{name}\n{}\n", String::from_utf8_lossy(bases)))
        .collect()
}

/// Records of a FASTA text with the bases of each joined on one line (gaps kept).
pub fn read_fasta(text: &str) -> Vec<Named> {
    let mut out: Vec<Named> = Vec::new();
    for line in text.lines() {
        if let Some(name) = line.strip_prefix('>') {
            out.push((name.trim().to_string(), Vec::new()));
        } else if let Some(last) = out.last_mut() {
            last.1.extend(line.trim().bytes());
        }
    }
    out
}

#[derive(Clone)]
pub struct Aligned {
    pub rows: Vec<Named>,
    pub engine: String,
    pub seconds: f64,
}

impl Aligned {
    pub fn columns(&self) -> usize {
        self.rows.first().map_or(0, |r| r.1.len())
    }

    pub fn msa(&self) -> Msa {
        Msa {
            names: self.rows.iter().map(|r| r.0.clone()).collect(),
            rows: self.rows.iter().map(|r| r.1.clone()).collect(),
        }
    }
}

/// Align `seqs`: with MAFFT when `mafft` is given, else with the built-in centre-star method.
pub fn align(seqs: &[Named], mafft: Option<&Installed>, dir: &Path) -> Result<Aligned, String> {
    if seqs.len() < 2 {
        return Err("Aligning needs at least two sequences.".into());
    }
    let start = Instant::now();
    let (rows, engine) = match mafft {
        Some(m) => {
            let named: Vec<Named> = ids(seqs.len())
                .into_iter()
                .zip(seqs.iter().map(|s| s.1.clone()))
                .collect();
            std::fs::write(dir.join("sequences.fasta"), fasta_text(&named))
                .map_err(|e| format!("Could not write the sequences: {e}"))?;
            let text = tools::run(&m.path, &["--auto", "--quiet", "sequences.fasta"], dir)?;
            let aligned = read_fasta(&text);
            if aligned.len() != seqs.len() {
                return Err(format!(
                    "MAFFT returned {} sequences for {}.",
                    aligned.len(),
                    seqs.len()
                ));
            }
            // MAFFT keeps the order; its names are the ids, and its letters come back in lower case
            let rows = aligned
                .into_iter()
                .zip(seqs)
                .map(|((_, bases), (name, _))| (name.clone(), bases.to_ascii_uppercase()))
                .collect();
            (rows, m.describe("MAFFT"))
        }
        None => {
            let msa = center_star(seqs)?;
            let rows = msa.names.into_iter().zip(msa.rows).collect();
            (rows, "centre-star alignment (built in)".to_string())
        }
    };
    Ok(Aligned {
        rows,
        engine,
        seconds: start.elapsed().as_secs_f64(),
    })
}

pub struct Tree {
    pub newick: String,
    pub engine: String,
    /// The substitution model chosen, as the program names it.
    pub model: String,
    pub ln_likelihood: Option<f64>,
    pub seconds: f64,
}

/// The text after `label` on the first line of `report` that has it.
fn after<'a>(report: &'a str, label: &str) -> Option<&'a str> {
    report
        .lines()
        .find_map(|l| l.split_once(label).map(|(_, rest)| rest.trim()))
}

/// A maximum-likelihood tree of the alignment: IQ-TREE (ModelFinder, and ultrafast bootstrap when
/// `bootstrap` is 1000 or more) when given, else the built-in search.
pub fn tree(
    aligned: &Aligned,
    iqtree: Option<&Installed>,
    bootstrap: usize,
    dir: &Path,
) -> Result<Tree, String> {
    let start = Instant::now();
    if aligned.rows.len() < 3 {
        return Err("A tree needs at least three sequences.".into());
    }
    let Some(program) = iqtree else {
        let found = crate::seq::likelihood::ml_tree(&aligned.msa())?;
        return Ok(Tree {
            newick: found.tree.newick(),
            engine: "likelihood search (built in)".into(),
            model: found.best.kind.label().to_string(),
            ln_likelihood: Some(found.best.ln_l),
            seconds: start.elapsed().as_secs_f64(),
        });
    };
    let named: Vec<Named> = ids(aligned.rows.len())
        .into_iter()
        .zip(aligned.rows.iter().map(|r| r.1.clone()))
        .collect();
    std::fs::write(dir.join("alignment.fasta"), fasta_text(&named))
        .map_err(|e| format!("Could not write the alignment: {e}"))?;
    let mut args = vec![
        "-s",
        "alignment.fasta",
        "-m",
        "MFP",
        "--prefix",
        "run",
        "-seed",
        "1",
        "-redo",
        "-quiet",
    ];
    let replicates = bootstrap.max(1000).to_string();
    if bootstrap > 0 {
        args.extend(["-B", &replicates]);
    }
    tools::run(&program.path, &args, dir)?;
    let treefile = std::fs::read_to_string(dir.join("run.treefile"))
        .map_err(|e| format!("IQ-TREE left no tree: {e}"))?;
    let report = std::fs::read_to_string(dir.join("run.iqtree")).unwrap_or_default();
    // the real names go back in place of s1, s2, ...
    let mut parsed = NwTree::parse_all(&treefile)?.remove(0);
    for node in &mut parsed.nodes {
        if let Some(k) = node
            .name
            .strip_prefix('s')
            .and_then(|n| n.parse::<usize>().ok())
        {
            if node.children.is_empty() {
                if let Some((name, _)) = aligned.rows.get(k - 1) {
                    node.name.clone_from(name);
                }
            }
        }
    }
    Ok(Tree {
        newick: parsed.newick(),
        engine: program.describe("IQ-TREE"),
        model: after(&report, "Best-fit model according to BIC:")
            .or_else(|| after(&report, "Model of substitution:"))
            .unwrap_or("")
            .to_string(),
        ln_likelihood: after(&report, "Log-likelihood of the tree:")
            .and_then(|t| t.split_whitespace().next()?.parse().ok()),
        seconds: start.elapsed().as_secs_f64(),
    })
}

/// The input primer3 reads (its "Boulder-IO" format) for `template`, with the program's settings as primer3
/// tags. `target` is the first base (from 0) and the length of the stretch the product must contain.
/// The dimer limits are the score-based ones (`PRIMER_MAX_SELF_ANY` and so on), which is what the built-in
/// rules apply; primer3 2.5 and newer would otherwise use its thermodynamic alignment, with other limits.
pub fn primer3_input(
    id: &str,
    template: &[u8],
    target: Option<(usize, usize)>,
    s: &Settings,
) -> String {
    let mut t = format!(
        "SEQUENCE_ID={id}\nSEQUENCE_TEMPLATE={}\n",
        String::from_utf8_lossy(template)
    );
    if let Some((first, len)) = target {
        t.push_str(&format!("SEQUENCE_TARGET={first},{len}\n"));
    }
    t.push_str(&format!(
        "PRIMER_TASK=generic\nPRIMER_PICK_LEFT_PRIMER=1\nPRIMER_PICK_RIGHT_PRIMER=1\n\
         PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0\nPRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=0\n\
         PRIMER_OPT_SIZE={}\nPRIMER_MIN_SIZE={}\nPRIMER_MAX_SIZE={}\n\
         PRIMER_OPT_TM={:.1}\nPRIMER_MIN_TM={:.1}\nPRIMER_MAX_TM={:.1}\n\
         PRIMER_MIN_GC={:.1}\nPRIMER_MAX_GC={:.1}\n\
         PRIMER_SALT_MONOVALENT={:.1}\nPRIMER_SALT_DIVALENT={:.1}\nPRIMER_DNTP_CONC={:.1}\nPRIMER_DNA_CONC={:.1}\n\
         PRIMER_MAX_POLY_X={}\nPRIMER_GC_CLAMP={}\nPRIMER_MAX_END_GC={}\n\
         PRIMER_MAX_SELF_ANY={}\nPRIMER_MAX_SELF_END={}\n\
         PRIMER_PAIR_MAX_COMPL_ANY={}\nPRIMER_PAIR_MAX_COMPL_END={}\nPRIMER_PAIR_MAX_DIFF_TM={:.1}\n\
         PRIMER_PRODUCT_SIZE_RANGE={}-{}\nPRIMER_NUM_RETURN={}\n=\n",
        s.opt_size, s.min_size, s.max_size,
        s.opt_tm, s.min_tm, s.max_tm,
        s.min_gc, s.max_gc,
        s.salt_mm, s.divalent_mm, s.dntp_mm, s.dna_nm,
        s.max_poly_x, s.gc_clamp, s.max_end_gc,
        s.max_self_any, s.max_self_end,
        s.pair_compl_any, s.pair_compl_end, s.max_diff_tm,
        s.product_min, s.product_max, s.num_return
    ));
    t
}

/// One pair as primer3 prints it.
#[derive(Clone, Debug, PartialEq)]
pub struct PrimerPair {
    pub left: String,
    pub right: String,
    pub left_tm: f64,
    pub right_tm: f64,
    /// The first template base each primer covers, from 0.
    pub left_start: usize,
    pub right_start: usize,
    pub product_size: usize,
    pub penalty: f64,
}

/// The pairs in primer3's output text.
pub fn read_primer3_output(text: &str) -> Vec<PrimerPair> {
    let get = |key: &str| {
        text.lines()
            .find_map(|l| l.strip_prefix(key)?.strip_prefix('='))
    };
    let number = |key: &str| get(key).and_then(|v| v.parse::<f64>().ok());
    (0..)
        .map_while(|i| {
            let left = get(&format!("PRIMER_LEFT_{i}_SEQUENCE"))?;
            let right = get(&format!("PRIMER_RIGHT_{i}_SEQUENCE"))?;
            let place = |key: &str| -> Option<(usize, usize)> {
                let (a, b) = get(key)?.split_once(',')?;
                Some((a.parse().ok()?, b.parse().ok()?))
            };
            let (left_start, _) = place(&format!("PRIMER_LEFT_{i}"))?;
            // primer3 gives the right primer's 5' end (its last template base) and its length
            let (right_last, right_len) = place(&format!("PRIMER_RIGHT_{i}"))?;
            Some(PrimerPair {
                left: left.to_string(),
                right: right.to_string(),
                left_tm: number(&format!("PRIMER_LEFT_{i}_TM"))?,
                right_tm: number(&format!("PRIMER_RIGHT_{i}_TM"))?,
                left_start,
                right_start: right_last + 1 - right_len,
                product_size: number(&format!("PRIMER_PAIR_{i}_PRODUCT_SIZE"))? as usize,
                penalty: number(&format!("PRIMER_PAIR_{i}_PENALTY"))?,
            })
        })
        .collect()
}

pub struct Primers {
    pub pairs: Vec<PrimerPair>,
    pub engine: String,
    pub seconds: f64,
}

/// Primer pairs for `template`: primer3 itself when given, else the built-in port of it.
pub fn primers(
    template: &[u8],
    target: Option<(usize, usize)>,
    settings: &Settings,
    primer3: Option<&Installed>,
    dir: &Path,
) -> Result<Primers, String> {
    let start = Instant::now();
    let (pairs, engine) = match primer3 {
        Some(p) => {
            let input = primer3_input("template", template, target, settings);
            std::fs::write(dir.join("primer3_input.txt"), input)
                .map_err(|e| format!("Could not write primer3's input: {e}"))?;
            let out = tools::run(&p.path, &["primer3_input.txt"], dir)?;
            (read_primer3_output(&out), p.describe("primer3"))
        }
        None => {
            let found = crate::seq::primer3::design(template, target, settings);
            let pairs = found
                .iter()
                .map(|p| PrimerPair {
                    left: String::from_utf8_lossy(&p.left.seq).into_owned(),
                    right: String::from_utf8_lossy(&p.right.seq).into_owned(),
                    left_tm: p.left.tm,
                    right_tm: p.right.tm,
                    left_start: p.left.start,
                    right_start: p.right.start,
                    product_size: p.product_size,
                    penalty: p.penalty,
                })
                .collect();
            (pairs, "primer3 rules (built in)".to_string())
        }
    };
    Ok(Primers {
        pairs,
        engine,
        seconds: start.elapsed().as_secs_f64(),
    })
}
