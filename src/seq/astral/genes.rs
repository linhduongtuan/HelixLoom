//! Gene trees for the quartet counting. Each tree is read, its taxa are numbered, and for every node that joins three or
//! more branches the taxa on each branch are kept as a set: that is all the counting needs to know of a gene tree.

use super::bits::Bits;
use crate::seq::newick::NwTree;
use std::collections::{BTreeSet, HashMap};

/// The fewest leaves a gene tree needs to say anything about quartets.
pub const MIN_LEAVES: usize = 4;

#[derive(Debug, Clone)]
pub struct GNode {
    pub parent: Option<usize>,
    pub children: Vec<usize>,
    pub taxon: Option<usize>,
}

/// A gene tree with its leaves numbered by taxon.
#[derive(Debug, Clone)]
pub struct GTree {
    pub nodes: Vec<GNode>,
    pub root: usize,
}

/// One branch at a node: the taxa on that side, and how many.
#[derive(Debug)]
pub struct Branch {
    pub set: Bits,
    pub size: u64,
}

/// What the counting reads of a gene tree: the taxa it has, and the branches around each node that joins at least three.
#[derive(Debug)]
pub struct Gene {
    pub present: Bits,
    pub nodes: Vec<Vec<Branch>>,
}

#[derive(Debug)]
pub struct Genes {
    /// The taxa, in name order; a taxon's number is its place here.
    pub taxa: Vec<String>,
    pub trees: Vec<GTree>,
    pub genes: Vec<Gene>,
    /// Trees left out for having fewer leaves than asked.
    pub skipped: usize,
}

impl GTree {
    fn from_newick(t: &NwTree, index: &HashMap<&str, usize>) -> GTree {
        let nodes = t
            .nodes
            .iter()
            .map(|n| GNode {
                parent: n.parent,
                children: n.children.clone(),
                taxon: n
                    .children
                    .is_empty()
                    .then(|| index.get(n.name.as_str()).copied())
                    .flatten(),
            })
            .collect();
        GTree {
            nodes,
            root: t.root,
        }
    }

    /// Children before parents.
    fn postorder(&self) -> Vec<usize> {
        let mut out = Vec::with_capacity(self.nodes.len());
        let mut stack = vec![(self.root, 0usize)];
        while let Some((i, next)) = stack.pop() {
            match self.nodes[i].children.get(next) {
                Some(&child) => {
                    stack.push((i, next + 1));
                    stack.push((child, 0));
                }
                None => out.push(i),
            }
        }
        out
    }

    /// The taxa below each node, by node number.
    pub fn clusters(&self, n: usize) -> Vec<Bits> {
        let mut sets = vec![Bits::new(n); self.nodes.len()];
        for i in self.postorder() {
            let node = &self.nodes[i];
            sets[i] = match node.taxon {
                Some(t) => Bits::of(n, [t]),
                None => node
                    .children
                    .iter()
                    .fold(Bits::new(n), |all, &c| all.union(&sets[c])),
            };
        }
        sets
    }

    /// Put a new leaf for `taxon` on the branch above `node` (above the root, a new root).
    pub fn insert_above(&mut self, node: usize, taxon: usize) {
        let leaf = self.nodes.len();
        let join = leaf + 1;
        let old_parent = self.nodes[node].parent;
        self.nodes.push(GNode {
            parent: Some(join),
            children: Vec::new(),
            taxon: Some(taxon),
        });
        self.nodes.push(GNode {
            parent: old_parent,
            children: vec![node, leaf],
            taxon: None,
        });
        self.nodes[node].parent = Some(join);
        match old_parent {
            Some(p) => {
                for c in &mut self.nodes[p].children {
                    if *c == node {
                        *c = join;
                    }
                }
            }
            None => self.root = join,
        }
    }

    pub fn gene(&self, n: usize) -> Gene {
        let sets = self.clusters(n);
        let present = sets[self.root].clone();
        let mut nodes = Vec::new();
        for (i, node) in self.nodes.iter().enumerate() {
            let mut branches: Vec<Branch> = node
                .children
                .iter()
                .map(|&c| &sets[c])
                .filter(|s| !s.is_empty())
                .map(|s| Branch {
                    set: s.clone(),
                    size: u64::from(s.count()),
                })
                .collect();
            let above = present.minus(&sets[i]);
            if !above.is_empty() {
                branches.push(Branch {
                    size: u64::from(above.count()),
                    set: above,
                });
            }
            if branches.len() >= 3 {
                nodes.push(branches);
            }
        }
        Gene { present, nodes }
    }
}

impl Gene {
    /// How many quartets of this tree it resolves: each is seen from both of the nodes that join its two pairs.
    pub fn resolved_quartets(&self) -> u128 {
        let mut twice: u128 = 0;
        for branches in &self.nodes {
            let total: u128 = branches.iter().map(|b| u128::from(b.size)).sum();
            let squares: u128 = branches.iter().map(|b| u128::from(b.size).pow(2)).sum();
            for b in branches {
                let k = u128::from(b.size);
                let pairs_in = k * k.saturating_sub(1) / 2;
                let rest = total - k;
                let across = (rest * rest - (squares - k * k)) / 2;
                twice += pairs_in * across;
            }
        }
        twice / 2
    }
}

fn leaf_names(t: &NwTree) -> Result<Vec<&str>, String> {
    let names: Vec<&str> = t
        .nodes
        .iter()
        .filter(|n| n.children.is_empty())
        .map(|n| n.name.as_str())
        .collect();
    if names.iter().any(|n| n.is_empty()) {
        return Err("A leaf of a gene tree has no name.".into());
    }
    Ok(names)
}

/// The most memory the sets of taxa of all the gene trees may take (see `set_bytes`): the search keeps them all.
const MAX_SET_BYTES: u128 = 8 << 30;

/// About what the sets of taxa of `nodes` tree nodes take for `taxa` taxa: a node keeps its own set and its children's, and
/// each set has a bit for every taxon.
fn set_bytes(nodes: u128, taxa: usize) -> u128 {
    nodes * 3 * (taxa.div_ceil(64) as u128) * 8
}

/// Whether the quartets of `genes` trees on `taxa` taxa can be counted in 64 bits: a tree has at most C(taxa, 4) of them
/// and every one is counted twice, once from each of the nodes that join its two pairs.
fn counts_fit(taxa: usize, genes: usize) -> bool {
    let n = taxa as u128;
    let quartets = n * n.saturating_sub(1) * n.saturating_sub(2) * n.saturating_sub(3) / 24;
    quartets * 2 * (genes as u128) < u128::from(u64::MAX)
}

/// Read gene trees (one `;`-ended Newick tree after the other). A tree with fewer than `min_leaves` leaves is left out;
/// a name that appears twice in one tree is refused (genes with several copies need ASTRAL-Pro, not built yet).
pub fn read(text: &str, min_leaves: usize) -> Result<Genes, String> {
    read_within(text, min_leaves, MAX_SET_BYTES)
}

/// `read`, refusing trees whose sets of taxa would take more than `max_bytes`, or whose quartets could not be counted.
fn read_within(text: &str, min_leaves: usize, max_bytes: u128) -> Result<Genes, String> {
    let parsed = NwTree::parse_all(text)?;
    let least = min_leaves.max(MIN_LEAVES);
    let mut kept = Vec::new();
    let mut names = BTreeSet::new();
    let mut skipped = 0;
    for (k, tree) in parsed.iter().enumerate() {
        let leaves = leaf_names(tree)?;
        let distinct: BTreeSet<&str> = leaves.iter().copied().collect();
        if distinct.len() != leaves.len() {
            return Err(format!(
                "Gene tree {} has the same leaf name more than once. Genes with several copies per species need ASTRAL-Pro, which is not built yet; give each leaf its species name only once per tree.",
                k + 1
            ));
        }
        if leaves.len() < least {
            skipped += 1;
            continue;
        }
        names.extend(leaves.iter().map(|s| s.to_string()));
        kept.push(tree);
    }
    if names.len() < MIN_LEAVES {
        return Err("A species tree needs at least four taxa in the gene trees.".into());
    }
    let taxa: Vec<String> = names.into_iter().collect();
    let nodes: u128 = kept.iter().map(|t| t.nodes.len() as u128).sum();
    let bytes = set_bytes(nodes, taxa.len());
    if bytes > max_bytes {
        let gigabytes = |b: u128| b.div_ceil(1 << 30);
        return Err(format!(
            "These {} gene trees on {} species are too big for the search: it would keep about {} GB of sets of species, and takes at most {} GB. Use fewer gene trees or fewer species.",
            kept.len(),
            taxa.len(),
            gigabytes(bytes),
            gigabytes(max_bytes)
        ));
    }
    if !counts_fit(taxa.len(), kept.len()) {
        return Err(format!(
            "There are too many species ({}) and gene trees ({}) for the search to count their quartets. Use fewer of them.",
            taxa.len(),
            kept.len()
        ));
    }
    let index: HashMap<&str, usize> = taxa
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    let trees: Vec<GTree> = kept.iter().map(|t| GTree::from_newick(t, &index)).collect();
    let genes = trees.iter().map(|t| t.gene(taxa.len())).collect();
    Ok(Genes {
        taxa,
        trees,
        genes,
        skipped,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_binary_tree_resolves_every_quartet_and_a_polytomy_fewer() {
        let g = read("((A,B),(C,(D,E)));\n(A,B,C,D,E);", 0).expect("trees");
        assert_eq!(g.taxa, ["A", "B", "C", "D", "E"]);
        assert_eq!(g.genes[0].resolved_quartets(), 5, "all C(5,4)");
        assert_eq!(g.genes[1].resolved_quartets(), 0, "a star resolves none");
        let partial = read("((A,B),C,D,E);", 0).expect("tree");
        // the cherry (A,B) with any two of C D E: three quartets
        assert_eq!(partial.genes[0].resolved_quartets(), 3);
    }

    #[test]
    fn a_rooted_binary_root_is_not_a_node_but_every_other_inner_node_is() {
        let g = read("((A,B),(C,D));", 0).expect("tree");
        // two inner nodes, each joining three branches (its two children and the other side)
        assert_eq!(g.genes[0].nodes.len(), 2);
        assert!(g.genes[0].nodes.iter().all(|b| b.len() == 3));
    }

    #[test]
    fn small_trees_are_skipped_and_repeated_names_refused() {
        let g = read("(A,B,C);\n((A,B),(C,D));", 0).expect("trees");
        assert_eq!((g.genes.len(), g.skipped), (1, 1));
        assert!(
            read("((A,B),(C,D));", 5).is_err(),
            "nothing is left when every tree is too small"
        );
        let err = read("((A,A),(C,D));", 0).expect_err("repeated name");
        assert!(err.contains("more than once"), "{err}");
    }

    #[test]
    fn trees_too_big_for_the_memory_or_for_the_counters_are_refused_with_what_to_do() {
        // ten trees of seven nodes on four taxa, each node with three sets of one word of eight bytes
        let text = "((A,B),(C,D));\n".repeat(10);
        let exactly = 10 * 7 * 3 * 8;
        assert!(read_within(&text, 0, exactly).is_ok());
        let why = read_within(&text, 0, exactly - 1).expect_err("too much memory");
        assert!(
            why.contains("fewer gene trees") && why.contains("10 gene trees on 4 species"),
            "{why}"
        );
        // C(10^4, 4) quartets, twice over, in a thousand trees fit in 64 bits; ten times the species do not
        assert!(counts_fit(10_000, 1_000));
        assert!(!counts_fit(100_000, 1_000));
        assert!(!counts_fit(30_000, 100_000));
        assert!(counts_fit(4, 1) && counts_fit(0, 0));
    }

    #[test]
    fn a_missing_taxon_can_be_put_back_on_a_branch() {
        let mut g = read("((A,B),(C,D));\n((A,B),(C,(D,E)));", 0).expect("trees");
        let n = g.taxa.len();
        assert_eq!(g.genes[0].present.count(), 4);
        let at = g.trees[0]
            .nodes
            .iter()
            .position(|x| x.taxon == Some(2))
            .expect("leaf C");
        g.trees[0].insert_above(at, 4);
        let again = g.trees[0].gene(n);
        assert_eq!(again.present.count(), 5);
        assert_eq!(again.resolved_quartets(), 5, "binary on five taxa");
    }
}
