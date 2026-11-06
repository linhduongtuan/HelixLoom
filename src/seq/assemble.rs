//! Contig assembler: overlap detection -> greedy layout -> consensus.
//!
//! Overlaps are ungapped (substitutions only). Reads are compared in both
//! orientations, joined best-overlap-first, and may close into a circle.
//! The overlap limits (min/max length, exact matching, circular product) follow the
//! original program's Gibson assembler defaults; the read-assembly scoring is our own.

use super::dna::{iupac_from_set, revcom};
use std::collections::HashMap;

/// Posting lists longer than this are skipped (low-complexity / repeat seeds).
const MAX_SEED_HITS: usize = 200;
/// Diagonals verified per (read, read, orientation) pair.
const MAX_DIAGONALS: usize = 3;
/// A mismatch costs this many matches when ranking overlaps.
const MISMATCH_COST: f64 = 3.0;

#[derive(Clone, Debug)]
pub struct Read {
    pub name: String,
    pub seq: Vec<u8>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Params {
    pub min_overlap: usize,
    pub min_identity: f64,
    pub max_overlap: Option<usize>,
    pub both_strands: bool,
    pub allow_containment: bool,
    pub circularize: bool,
    pub only_circular: bool,
}

impl Default for Params {
    fn default() -> Self {
        Self::reads()
    }
}

impl Params {
    /// Sequencing reads: tolerant overlaps, linear contigs.
    pub fn reads() -> Self {
        Params {
            min_overlap: 20,
            min_identity: 0.90,
            max_overlap: None,
            both_strands: true,
            allow_containment: true,
            circularize: true,
            only_circular: false,
        }
    }

    /// Gibson-style fragments: exact 12-100 bp end overlaps, circular product only.
    pub fn gibson() -> Self {
        Params {
            min_overlap: 12,
            min_identity: 1.0,
            max_overlap: Some(100),
            both_strands: true,
            allow_containment: false,
            circularize: true,
            only_circular: true,
        }
    }

    fn kmer(&self) -> usize {
        self.min_overlap.clamp(4, 12)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Placed {
    pub read: usize,
    pub reverse: bool,
    pub offset: usize,
    pub contained: bool,
}

#[derive(Clone, Debug)]
pub struct Contig {
    pub name: String,
    pub placed: Vec<Placed>,
    pub consensus: Vec<u8>,
    /// Columns covered by reads; equals `consensus.len()` unless circular.
    pub span: usize,
    pub circular: bool,
    pub closing_overlap: usize,
    pub depth: Vec<u32>,
    pub conflicts: usize,
    /// Consensus columns where at least one read disagrees with the majority.
    pub disagree: Vec<bool>,
}

#[derive(Clone, Debug, Default)]
pub struct Assembly {
    pub contigs: Vec<Contig>,
    pub unassembled: Vec<usize>,
    pub n50: usize,
    pub total_length: usize,
}

impl Contig {
    pub fn mean_depth(&self) -> f64 {
        if self.depth.is_empty() {
            return 0.0;
        }
        self.depth.iter().map(|&d| d as f64).sum::<f64>() / self.depth.len() as f64
    }

    /// Consensus base at a layout column (circular contigs wrap).
    pub fn consensus_at(&self, col: usize) -> u8 {
        self.consensus[col % self.consensus.len().max(1)]
    }
}

// ---------------------------------------------------------------------------
// Overlap detection
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
struct Ov {
    a: usize,
    b: usize,
    b_rev: bool,
    /// `b` starts at `off` in `a`'s coordinates (negative: starts before `a`).
    off: isize,
    len: usize,
    score: f64,
}

#[derive(Clone, Copy)]
struct Rel {
    container: usize,
    off: isize,
    member_rev: bool,
    score: f64,
}

#[derive(Clone, Copy)]
struct Edge {
    from: usize,
    to: usize,
    off: usize,
    len: usize,
    score: f64,
}

/// Upper-case A/C/G/T with everything else as N (gaps and whitespace dropped, U read as T).
pub fn normalize(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .filter(|c| !matches!(c, b'-' | b'.' | b' ' | b'\n' | b'\r' | b'\t'))
        .map(|&c| match c.to_ascii_uppercase() {
            b'U' => b'T',
            b @ (b'A' | b'C' | b'G' | b'T') => b,
            _ => b'N',
        })
        .collect()
}

fn base_code(c: u8) -> Option<u32> {
    match c {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Count matching positions of `b` placed at `d` on `a`. `None` when they don't overlap.
fn verify(a: &[u8], b: &[u8], d: isize) -> Option<(usize, usize)> {
    let start = d.max(0) as usize;
    let end = (a.len() as isize).min(d + b.len() as isize);
    if end <= start as isize {
        return None;
    }
    let end = end as usize;
    let matches = (start..end)
        .filter(|&x| {
            let (ca, cb) = (a[x], b[(x as isize - d) as usize]);
            ca == cb || ca == b'N' || cb == b'N'
        })
        .count();
    Some((end - start, matches))
}

/// Overlaps whose length/identity pass `p`; `None` when the placement isn't an end overlap or containment.
fn judge(a: &[u8], b: &[u8], d: isize, p: &Params) -> Option<(usize, f64)> {
    let (la, lb) = (a.len() as isize, b.len() as isize);
    let contained = (d >= 0 && d + lb <= la) || (d <= 0 && d + lb >= la);
    if contained && !p.allow_containment {
        return None;
    }
    if !contained && p.max_overlap.is_some_and(|m| overlap_len(la, lb, d) > m) {
        return None;
    }
    let (len, matches) = verify(a, b, d)?;
    if len < p.min_overlap || (matches as f64) < p.min_identity * len as f64 - 1e-9 {
        return None;
    }
    let mism = (len - matches) as f64;
    Some((len, matches as f64 - MISMATCH_COST * mism))
}

fn overlap_len(la: isize, lb: isize, d: isize) -> usize {
    (la.min(d + lb) - d.max(0)).max(0) as usize
}

fn find_overlaps(fwd: &[Vec<u8>], rev: &[Vec<u8>], p: &Params) -> Vec<Ov> {
    let k = p.kmer();
    let mut index: HashMap<u32, Vec<(u32, u32)>> = HashMap::new();
    for (i, s) in fwd.iter().enumerate() {
        for (pos, code) in kmers(s, k) {
            index.entry(code).or_default().push((i as u32, pos as u32));
        }
    }
    let mut found = Vec::new();
    for j in 0..fwd.len() {
        let orientations: &[bool] = if p.both_strands {
            &[false, true]
        } else {
            &[false]
        };
        for &b_rev in orientations {
            let b = if b_rev { &rev[j] } else { &fwd[j] };
            let mut votes: HashMap<(usize, isize), u32> = HashMap::new();
            for (pos, code) in kmers(b, k) {
                let Some(hits) = index.get(&code) else {
                    continue;
                };
                if hits.len() > MAX_SEED_HITS {
                    continue;
                }
                for &(i, pi) in hits.iter().filter(|(i, _)| (*i as usize) < j) {
                    *votes
                        .entry((i as usize, pi as isize - pos as isize))
                        .or_default() += 1;
                }
            }
            collect_best(&mut found, votes, fwd, b, j, b_rev, p);
        }
    }
    found
}

fn collect_best(
    found: &mut Vec<Ov>,
    votes: HashMap<(usize, isize), u32>,
    fwd: &[Vec<u8>],
    b: &[u8],
    j: usize,
    b_rev: bool,
    p: &Params,
) {
    let mut by_read: HashMap<usize, Vec<(isize, u32)>> = HashMap::new();
    for ((i, d), n) in votes {
        by_read.entry(i).or_default().push((d, n));
    }
    for (i, mut diags) in by_read {
        diags.sort_by(|x, y| y.1.cmp(&x.1).then(x.0.cmp(&y.0)));
        let best = diags
            .iter()
            .take(MAX_DIAGONALS)
            .filter_map(|&(d, _)| judge(&fwd[i], b, d, p).map(|(len, score)| (d, len, score)))
            .max_by(|x, y| x.2.total_cmp(&y.2).then(y.0.cmp(&x.0)));
        if let Some((off, len, score)) = best {
            found.push(Ov {
                a: i,
                b: j,
                b_rev,
                off,
                len,
                score,
            });
        }
    }
}

fn kmers(seq: &[u8], k: usize) -> impl Iterator<Item = (usize, u32)> + '_ {
    let mask = if k >= 16 {
        u32::MAX
    } else {
        (1u32 << (2 * k)) - 1
    };
    let mut code = 0u32;
    let mut valid = 0usize;
    seq.iter().enumerate().filter_map(move |(pos, &c)| {
        match base_code(c) {
            Some(v) => {
                code = ((code << 2) | v) & mask;
                valid += 1;
            }
            None => valid = 0,
        }
        (valid >= k).then(|| (pos + 1 - k, code))
    })
}

// ---------------------------------------------------------------------------
// Layout
// ---------------------------------------------------------------------------

fn id(read: usize, rev: bool) -> usize {
    read * 2 + rev as usize
}

/// Oriented-read graph: `out[x] = (next, offset, overlap)`, `inn[x] = previous`.
struct Layout {
    out: Vec<Option<(usize, usize, usize)>>,
    inn: Vec<Option<usize>>,
}

fn mirror(e: Edge, lens: &[usize]) -> Edge {
    let (lx, ly) = (lens[e.from / 2], lens[e.to / 2]);
    Edge {
        from: e.to ^ 1,
        to: e.from ^ 1,
        off: e.off + ly - lx,
        len: e.len,
        score: e.score,
    }
}

fn dovetail_edges(ovs: &[Ov], contained: &[bool], lens: &[usize]) -> Vec<Edge> {
    let mut edges = Vec::new();
    for o in ovs.iter().filter(|o| !contained[o.a] && !contained[o.b]) {
        let (la, lb) = (lens[o.a] as isize, lens[o.b] as isize);
        let is_containment = (o.off >= 0 && o.off + lb <= la) || (o.off <= 0 && o.off + lb >= la);
        if is_containment {
            continue;
        }
        let (a, b) = (id(o.a, false), id(o.b, o.b_rev));
        let e = if o.off > 0 {
            Edge {
                from: a,
                to: b,
                off: o.off as usize,
                len: o.len,
                score: o.score,
            }
        } else {
            Edge {
                from: b,
                to: a,
                off: (-o.off) as usize,
                len: o.len,
                score: o.score,
            }
        };
        edges.push(e);
        edges.push(mirror(e, lens));
    }
    edges.sort_by(|x, y| {
        y.score
            .total_cmp(&x.score)
            .then(y.len.cmp(&x.len))
            .then(x.from.cmp(&y.from))
            .then(x.to.cmp(&y.to))
    });
    edges
}

fn greedy_layout(lens: &[usize], edges: &[Edge], circularize: bool) -> Layout {
    let n = lens.len();
    let mut lay = Layout {
        out: vec![None; 2 * n],
        inn: vec![None; 2 * n],
    };
    let mut uf: Vec<usize> = (0..n).collect();
    for e in edges {
        if lay.out[e.from].is_some() || lay.inn[e.to].is_some() || e.from / 2 == e.to / 2 {
            continue;
        }
        let (ra, rb) = (find(&mut uf, e.from / 2), find(&mut uf, e.to / 2));
        if ra == rb && !circularize {
            continue;
        }
        for x in [*e, mirror(*e, lens)] {
            lay.out[x.from] = Some((x.to, x.off, x.len));
            lay.inn[x.to] = Some(x.from);
        }
        if ra != rb {
            uf[ra] = rb;
        }
    }
    lay
}

fn find(uf: &mut [usize], mut x: usize) -> usize {
    while uf[x] != x {
        uf[x] = uf[uf[x]];
        x = uf[x];
    }
    x
}

/// Containment relations to top-level (not themselves contained) reads, best score first.
fn containment(ovs: &[Ov], lens: &[usize]) -> Vec<Option<Rel>> {
    let mut rels: Vec<Vec<Rel>> = vec![Vec::new(); lens.len()];
    for o in ovs {
        let (la, lb) = (lens[o.a] as isize, lens[o.b] as isize);
        if o.off >= 0 && o.off + lb <= la {
            rels[o.b].push(Rel {
                container: o.a,
                off: o.off,
                member_rev: o.b_rev,
                score: o.score,
            });
        } else if o.off <= 0 && o.off + lb >= la {
            let (off, member_rev) = if o.b_rev {
                (lb + o.off - la, true)
            } else {
                (-o.off, false)
            };
            rels[o.a].push(Rel {
                container: o.b,
                off,
                member_rev,
                score: o.score,
            });
        }
    }
    // A read may only be contained by a strictly longer read, or an equal one with a lower index.
    let wins = |c: usize, m: usize| lens[c] > lens[m] || (lens[c] == lens[m] && c < m);
    let mut order: Vec<usize> = (0..lens.len()).collect();
    order.sort_by(|&x, &y| lens[y].cmp(&lens[x]).then(x.cmp(&y)));
    let mut chosen: Vec<Option<Rel>> = vec![None; lens.len()];
    for &m in &order {
        chosen[m] = rels[m]
            .iter()
            .filter(|r| wins(r.container, m) && chosen[r.container].is_none())
            .max_by(|x, y| {
                x.score
                    .total_cmp(&y.score)
                    .then(y.container.cmp(&x.container))
            })
            .copied();
    }
    chosen
}

// ---------------------------------------------------------------------------
// Assembly driver
// ---------------------------------------------------------------------------

pub fn assemble(reads: &[Read], p: &Params) -> Assembly {
    let fwd: Vec<Vec<u8>> = reads.iter().map(|r| normalize(&r.seq)).collect();
    let rev: Vec<Vec<u8>> = fwd.iter().map(|s| revcom(s)).collect();
    let lens: Vec<usize> = fwd.iter().map(Vec::len).collect();

    let ovs: Vec<Ov> = if reads.len() > 1 {
        find_overlaps(&fwd, &rev, p)
    } else {
        Vec::new()
    };
    let rels = containment(&ovs, &lens);
    let contained: Vec<bool> = rels.iter().map(Option::is_some).collect();

    let edges = dovetail_edges(&ovs, &contained, &lens);
    let lay = greedy_layout(&lens, &edges, p.circularize);

    let mut asm = Assembly::default();
    let mut visited = vec![false; reads.len()];
    let mut members: HashMap<usize, Vec<usize>> = HashMap::new();
    for (m, r) in rels.iter().enumerate() {
        if let Some(r) = r {
            members.entry(r.container).or_default().push(m);
        }
    }
    for r in 0..reads.len() {
        if visited[r] || contained[r] {
            continue;
        }
        let chain = walk_chain(&lay, r, &mut visited);
        let mut placed: Vec<Placed> = chain
            .nodes
            .iter()
            .map(|&(rid, off)| Placed {
                read: rid / 2,
                reverse: rid % 2 == 1,
                offset: off,
                contained: false,
            })
            .collect();
        for node in placed.clone() {
            for &m in members.get(&node.read).into_iter().flatten() {
                if let Some(rel) = &rels[m] {
                    placed.push(place_member(rel, m, &node, &lens));
                }
            }
        }
        // A lone read only becomes a contig when its own ends overlap into a circle.
        let lone_circle = (placed.len() == 1 && p.circularize)
            .then(|| self_overlap(&fwd[r], p))
            .flatten()
            .map(|ov| (lens[r] - ov, ov));
        if placed.len() == 1 && lone_circle.is_none() {
            asm.unassembled.push(r);
            continue;
        }
        asm.contigs.push(build_contig(
            placed,
            chain.circular_close.or(lone_circle),
            &fwd,
            &rev,
        ));
    }
    finish(asm, p)
}

struct Chain {
    nodes: Vec<(usize, usize)>,
    circular_close: Option<(usize, usize)>,
}

fn walk_chain(lay: &Layout, start: usize, visited: &mut [bool]) -> Chain {
    let mut head = id(start, false);
    let mut steps = 0;
    while let Some(prev) = lay.inn[head] {
        head = prev;
        steps += 1;
        if head == id(start, false) || steps > lay.inn.len() {
            break;
        }
    }
    let mut nodes = vec![(head, 0usize)];
    visited[head / 2] = true;
    let mut cur = head;
    let mut off = 0usize;
    let mut close = None;
    while let Some((next, o, len)) = lay.out[cur] {
        if next == head {
            close = Some((off + o, len));
            break;
        }
        if visited[next / 2] {
            break;
        }
        off += o;
        nodes.push((next, off));
        visited[next / 2] = true;
        cur = next;
    }
    Chain {
        nodes,
        circular_close: close,
    }
}

fn place_member(rel: &Rel, member: usize, container: &Placed, lens: &[usize]) -> Placed {
    let lc = lens[container.read] as isize;
    let lm = lens[member] as isize;
    let (pos, member_rev) = if container.reverse {
        (
            container.offset as isize + (lc - (rel.off + lm)),
            !rel.member_rev,
        )
    } else {
        (container.offset as isize + rel.off, rel.member_rev)
    };
    Placed {
        read: member,
        reverse: member_rev,
        offset: pos.max(0) as usize,
        contained: true,
    }
}

/// Longest end overlap of a read with itself (circular closure of one fragment).
fn self_overlap(seq: &[u8], p: &Params) -> Option<usize> {
    let len = seq.len();
    (p.min_overlap..=len / 2).rev().find(|&l| {
        if p.max_overlap.is_some_and(|m| l > m) {
            return false;
        }
        let m = (0..l)
            .filter(|&x| seq[x] == seq[len - l + x] || seq[x] == b'N' || seq[len - l + x] == b'N')
            .count();
        m as f64 >= p.min_identity * l as f64 - 1e-9
    })
}

fn build_contig(
    mut placed: Vec<Placed>,
    closing: Option<(usize, usize)>,
    fwd: &[Vec<u8>],
    rev: &[Vec<u8>],
) -> Contig {
    let min_off = placed.iter().map(|p| p.offset).min().unwrap_or(0);
    placed.iter_mut().for_each(|p| p.offset -= min_off);
    placed.sort_by_key(|p| (p.offset, p.contained, p.read));
    let seq_of = |p: &Placed| {
        if p.reverse {
            &rev[p.read]
        } else {
            &fwd[p.read]
        }
    };
    let span = placed
        .iter()
        .map(|p| p.offset + seq_of(p).len())
        .max()
        .unwrap_or(0);
    let circ_len = closing.map(|(c, _)| c).filter(|&c| c > 0 && c < span);
    let width = circ_len.unwrap_or(span);

    let mut counts = vec![[0u32; 4]; width];
    let mut depth = vec![0u32; span];
    for p in &placed {
        for (k, &c) in seq_of(p).iter().enumerate() {
            let x = p.offset + k;
            depth[x] += 1;
            if let Some(v) = base_code(c) {
                counts[x % width][v as usize] += 1;
            }
        }
    }
    let (consensus, disagree) = call_consensus(&counts);
    Contig {
        name: String::new(),
        placed,
        consensus,
        span,
        circular: circ_len.is_some(),
        closing_overlap: closing.filter(|_| circ_len.is_some()).map_or(0, |(_, l)| l),
        depth,
        conflicts: disagree.iter().filter(|&&d| d).count(),
        disagree,
    }
}

/// Majority base per column (ties become an IUPAC code) and whether any read disagreed.
fn call_consensus(counts: &[[u32; 4]]) -> (Vec<u8>, Vec<bool>) {
    counts
        .iter()
        .map(|c| {
            let top = *c.iter().max().unwrap_or(&0);
            let total: u32 = c.iter().sum();
            if total == 0 {
                return (b'N', false);
            }
            let tied: Vec<u8> = (0..4)
                .filter(|&i| c[i] == top)
                .map(|i| b"ACGT"[i])
                .collect();
            let base = if tied.len() == 1 {
                tied[0]
            } else {
                iupac_from_set(&tied)
            };
            (base, top < total)
        })
        .unzip()
}

fn finish(mut asm: Assembly, p: &Params) -> Assembly {
    if p.only_circular {
        let dropped: Vec<usize> = asm
            .contigs
            .iter()
            .filter(|c| !c.circular)
            .flat_map(|c| c.placed.iter().map(|x| x.read))
            .collect();
        asm.contigs.retain(|c| c.circular);
        asm.unassembled.extend(dropped);
    }
    asm.contigs.sort_by(|a, b| {
        b.consensus
            .len()
            .cmp(&a.consensus.len())
            .then(a.placed[0].read.cmp(&b.placed[0].read))
    });
    for (i, c) in asm.contigs.iter_mut().enumerate() {
        c.name = format!("Contig{}", i + 1);
    }
    asm.unassembled.sort_unstable();
    asm.unassembled.dedup();
    asm.total_length = asm.contigs.iter().map(|c| c.consensus.len()).sum();
    asm.n50 = n50(
        asm.contigs.iter().map(|c| c.consensus.len()).collect(),
        asm.total_length,
    );
    asm
}

fn n50(mut sizes: Vec<usize>, total: usize) -> usize {
    sizes.sort_unstable_by(|a, b| b.cmp(a));
    let mut acc = 0;
    for s in sizes {
        acc += s;
        if acc * 2 >= total {
            return s;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg_genome(len: usize, mut seed: u64) -> Vec<u8> {
        (0..len)
            .map(|_| {
                seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[((seed >> 33) % 4) as usize]
            })
            .collect()
    }

    fn reads_from(genome: &[u8], width: usize, step: usize, flip_every: usize) -> Vec<Read> {
        let mut out = Vec::new();
        let mut start = 0;
        while start < genome.len() {
            let end = (start + width).min(genome.len());
            let mut s = genome[start..end].to_vec();
            if out.len() % flip_every == 1 {
                s = revcom(&s);
            }
            out.push(Read {
                name: format!("r{}", out.len()),
                seq: s,
            });
            if end == genome.len() {
                break;
            }
            start += step;
        }
        out
    }

    fn same_or_revcom(a: &[u8], b: &[u8]) -> bool {
        a == b || revcom(a) == b
    }

    #[test]
    fn assembles_mixed_strand_reads_into_one_contig() {
        let g = lcg_genome(3000, 7);
        let reads = reads_from(&g, 400, 250, 2);
        let asm = assemble(&reads, &Params::reads());
        assert_eq!(asm.contigs.len(), 1);
        assert!(asm.unassembled.is_empty());
        assert!(same_or_revcom(&asm.contigs[0].consensus, &g));
        assert!(!asm.contigs[0].circular);
    }

    #[test]
    fn votes_away_sequencing_errors() {
        let g = lcg_genome(2400, 11);
        let mut reads = reads_from(&g, 500, 100, 3);
        let count = reads.len();
        // One substitution per interior read, at a genome position covered by 3+ reads.
        for (n, r) in reads
            .iter_mut()
            .enumerate()
            .filter(|(n, _)| (2..count - 2).contains(n))
        {
            let gp = 100 * n + 150 + (n * 7) % 200;
            let idx = if n % 3 == 1 {
                r.seq.len() - 1 - (gp - 100 * n)
            } else {
                gp - 100 * n
            };
            r.seq[idx] = if r.seq[idx] == b'A' { b'C' } else { b'A' };
        }
        let asm = assemble(&reads, &Params::reads());
        assert_eq!(asm.contigs.len(), 1);
        assert!(same_or_revcom(&asm.contigs[0].consensus, &g));
        assert!(asm.contigs[0].conflicts > 0);
        assert!(asm.contigs[0].mean_depth() > 3.0);
    }

    #[test]
    fn closes_overlapping_fragments_into_a_circle() {
        let g = lcg_genome(1000, 3);
        let frag = |a: usize, b: usize| {
            g.iter()
                .cycle()
                .skip(a)
                .take(b - a)
                .copied()
                .collect::<Vec<u8>>()
        };
        let reads: Vec<Read> = [(0, 360), (330, 690), (660, 1030)]
            .iter()
            .enumerate()
            .map(|(i, &(a, b))| Read {
                name: format!("f{i}"),
                seq: frag(a, b),
            })
            .collect();
        let asm = assemble(&reads, &Params::gibson());
        assert_eq!(asm.contigs.len(), 1);
        let c = &asm.contigs[0];
        assert!(c.circular);
        assert_eq!(c.consensus.len(), 1000);
        let doubled: Vec<u8> = g.iter().chain(g.iter()).copied().collect();
        let rc_doubled = revcom(&doubled);
        let found = |hay: &[u8]| hay.windows(1000).any(|w| w == c.consensus.as_slice());
        assert!(found(&doubled) || found(&rc_doubled));
    }

    #[test]
    fn places_contained_reads_inside_their_container() {
        let g = lcg_genome(900, 5);
        let mut reads = reads_from(&g, 500, 400, 5);
        reads.push(Read {
            name: "inner".into(),
            seq: revcom(&g[100..220]),
        });
        let asm = assemble(&reads, &Params::reads());
        assert_eq!(asm.contigs.len(), 1);
        let inner = asm.contigs[0]
            .placed
            .iter()
            .find(|p| p.read == reads.len() - 1)
            .unwrap();
        assert!(inner.contained);
        assert!(same_or_revcom(&asm.contigs[0].consensus, &g));
    }

    #[test]
    fn keeps_unrelated_reads_apart() {
        let a = lcg_genome(1200, 21);
        let b = lcg_genome(1200, 99);
        let mut reads = reads_from(&a, 400, 300, 2);
        reads.extend(reads_from(&b, 400, 300, 2));
        reads.push(Read {
            name: "stray".into(),
            seq: lcg_genome(300, 1234),
        });
        let asm = assemble(&reads, &Params::reads());
        assert_eq!(asm.contigs.len(), 2);
        assert_eq!(asm.unassembled, vec![reads.len() - 1]);
        assert_eq!(asm.total_length, 2400);
        assert_eq!(asm.n50, 1200);
    }

    #[test]
    fn exact_mode_rejects_overlaps_with_mismatches() {
        let g = lcg_genome(700, 8);
        let mut reads = reads_from(&g, 400, 300, 5);
        // The second read is stored reverse-complemented, so its overlap with the first is at its end.
        let at = reads[1].seq.len() - 50;
        reads[1].seq[at] = if reads[1].seq[at] == b'G' { b'T' } else { b'G' };
        let strict = Params {
            min_identity: 1.0,
            ..Params::reads()
        };
        assert_eq!(assemble(&reads, &strict).contigs.len(), 0);
        assert_eq!(assemble(&reads, &Params::reads()).contigs.len(), 1);
    }
}
