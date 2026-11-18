//! What one evaluation of the likelihood of a tree needs of the tree, worked out once on the host, whichever library does the tensor
//! work (`gpu`, candle; `torch`, PyTorch): the order in which the inner nodes are done and the transition matrix of every branch.

use super::{expm, Matrix, Tree};

/// One side of an inner node: a tip (its number) or an inner node (its number), and the place of the matrix of its branch.
pub(super) struct Child {
    pub tip: Option<usize>,
    pub node: usize,
    pub slot: usize,
}

pub(super) struct Step {
    pub node: usize,
    pub children: [Child; 2],
}

/// What one evaluation needs of the tree, worked out once on the host: the order in which the inner nodes are done, and the
/// transition matrix of every branch (in 32 bits: sixteen numbers for a branch to an inner node, a table of five rows of four
/// for a branch to a tip).
pub(super) struct Plan {
    pub steps: Vec<Step>,
    pub first: usize,
    pub inner_p: Vec<f32>,
    pub tip_p: Vec<f32>,
    pub root_p: [f32; 16],
    /// Five rows of four: the base frequencies at the place of each base (and all of them for "any base"), what tip 0 says at the root.
    pub root_table: [f32; 20],
    pub nodes: usize,
}

/// The steps of a pruning, the matrices of the branches to inner nodes, and the tables of the branches to tips.
type Pruned = (Vec<Step>, Vec<f32>, Vec<f32>);

/// The order in which the inner nodes of the subtree at `first`, away from its neighbour `parent`, are done, and the transition matrices of
/// its branches (in 32 bits, as in `Plan`).
fn prune_steps(tree: &Tree, first: usize, parent: usize, q: &Matrix) -> Result<Pruned, String> {
    // The inner nodes below `first`, children before parents: go down with a stack and read the order backwards (the child pushed
    // first is the one done first). A node's vector is kept until its parent is done, so doing the larger subtree first, which needs
    // the more vectors, holds one vector while the smaller is done: at most about the logarithm of the tips' number of them at once,
    // where the other order can hold one for every level of a tree that is a spine with a cherry at each level (the card's memory).
    let children_of = |node: usize, parent: usize| -> Vec<usize> {
        tree.adj[node]
            .iter()
            .filter(|e| e.0 != parent && e.0 >= tree.tips)
            .map(|e| e.0)
            .collect()
    };
    let mut visit: Vec<(usize, usize)> = Vec::new();
    let mut down = vec![(first, parent)];
    while let Some((node, parent)) = down.pop() {
        visit.push((node, parent));
        down.extend(children_of(node, parent).into_iter().map(|c| (c, node)));
    }
    // how many vectors each subtree needs: the larger need first, then one more for the vector held meanwhile
    let mut need = vec![0usize; tree.adj.len()];
    for &(node, parent) in visit.iter().rev() {
        let mut inner: Vec<usize> = children_of(node, parent).iter().map(|&c| need[c]).collect();
        inner.sort_unstable_by(|a, b| b.cmp(a));
        need[node] = match inner.as_slice() {
            [] => 1,
            [one] => (*one).max(1),
            [a, b, ..] => (*a).max(*b + 1),
        };
    }
    let mut down = vec![(first, parent)];
    let mut reverse: Vec<(usize, usize)> = Vec::new();
    while let Some((node, parent)) = down.pop() {
        reverse.push((node, parent));
        let mut inner = children_of(node, parent);
        inner.sort_by(|a, b| need[*b].cmp(&need[*a]));
        down.extend(inner.into_iter().map(|c| (c, node)));
    }
    let (mut inner_p, mut tip_p) = (Vec::new(), Vec::new());
    let mut steps = Vec::with_capacity(reverse.len());
    for &(node, parent) in reverse.iter().rev() {
        let children: Vec<&(usize, f64)> =
            tree.adj[node].iter().filter(|e| e.0 != parent).collect();
        let [a, b] = children.as_slice() else {
            return Err("An inner node of the tree does not have two children.".into());
        };
        let mut make = |&&(c, length): &&(usize, f64)| -> Child {
            let p = expm(q, length);
            if c < tree.tips {
                let slot = tip_p.len() / 20;
                // the table: for a base, the chance of the parent's every base given that one at the tip (a column of P)
                for s in 0..4 {
                    tip_p.extend(p.iter().map(|row| row[s] as f32));
                }
                tip_p.extend([1.0f32; 4]);
                Child {
                    tip: Some(c),
                    node: c,
                    slot,
                }
            } else {
                let slot = inner_p.len() / 16;
                inner_p.extend(p.iter().flatten().map(|&x| x as f32));
                Child {
                    tip: None,
                    node: c,
                    slot,
                }
            }
        };
        steps.push(Step {
            node,
            children: [make(a), make(b)],
        });
    }
    Ok((steps, inner_p, tip_p))
}

/// What a worker needs to give, for each pattern, the likelihood of what lies beyond one end of a branch given each state at that end:
/// the pruning of the subtree at `node` away from its neighbour `from`. A tip is a lookup in a table and needs no steps.
pub(super) struct Beyond {
    pub steps: Vec<Step>,
    /// The tip, or the inner node whose vector is the answer.
    pub node: usize,
    pub is_tip: bool,
    pub inner_p: Vec<f32>,
    pub tip_p: Vec<f32>,
}

impl Beyond {
    pub(super) fn new(tree: &Tree, node: usize, from: usize, q: &Matrix) -> Result<Beyond, String> {
        if node < tree.tips {
            return Ok(Beyond {
                steps: Vec::new(),
                node,
                is_tip: true,
                inner_p: Vec::new(),
                tip_p: Vec::new(),
            });
        }
        let (steps, inner_p, tip_p) = prune_steps(tree, node, from, q)?;
        Ok(Beyond {
            steps,
            node,
            is_tip: false,
            inner_p,
            tip_p,
        })
    }

    pub(super) fn steps_table(&self) -> Vec<u32> {
        steps_table(&self.steps)
    }
}

/// The steps as numbers for a worker in another language, seven for each: the node, then for each of its two children whether it is a
/// tip (1) or an inner node (0), its number, and the place of its matrix.
fn steps_table(steps: &[Step]) -> Vec<u32> {
    steps
        .iter()
        .flat_map(|step| {
            let [a, b] = &step.children;
            [
                step.node as u32,
                u32::from(a.tip.is_some()),
                a.node as u32,
                a.slot as u32,
                u32::from(b.tip.is_some()),
                b.node as u32,
                b.slot as u32,
            ]
        })
        .collect()
}

impl Plan {
    pub(super) fn new(tree: &Tree, q: &Matrix, freqs: &[f64; 4]) -> Result<Plan, String> {
        let (first, root_length) = *tree
            .adj
            .first()
            .and_then(|a| a.first())
            .ok_or("The tree has no branch.")?;
        if first < tree.tips {
            return Err("The tree has fewer than three tips.".into());
        }
        let (steps, inner_p, tip_p) = prune_steps(tree, first, 0, q)?;
        let p = expm(q, root_length);
        let mut root_p = [0.0f32; 16];
        for (slot, x) in root_p.iter_mut().zip(p.iter().flatten()) {
            *slot = *x as f32;
        }
        let mut root_table = [0.0f32; 20];
        for s in 0..4 {
            root_table[s * 4 + s] = freqs[s] as f32;
        }
        for i in 0..4 {
            root_table[16 + i] = freqs[i] as f32;
        }
        Ok(Plan {
            steps,
            first,
            inner_p,
            tip_p,
            root_p,
            root_table,
            nodes: tree.adj.len(),
        })
    }

    /// The steps as numbers for a worker in another language (see `steps_table`).
    pub(super) fn steps_table(&self) -> Vec<u32> {
        steps_table(&self.steps)
    }
}

#[cfg(test)]
mod tests {
    use super::super::rate_matrix;
    use super::*;

    /// The most vectors that are held at once when the steps are done in order: a node's vector is kept until its parent is done.
    fn peak_live(plan: &Plan) -> usize {
        let (mut live, mut peak) = (0usize, 0usize);
        for step in &plan.steps {
            peak = peak.max(live);
            live -= step.children.iter().filter(|c| c.tip.is_none()).count();
            live += 1;
            peak = peak.max(live);
        }
        peak
    }

    /// Tip 0, then a spine of `cherries - 1` nodes, each with a cherry (an inner node with two tips) and the rest of the spine as its
    /// children, the cherry listed first; the last spine node has two cherries. A tree that kept a cherry's vector while it did the
    /// whole spine would hold one for every level.
    fn spine_of_cherries(cherries: usize, cherry_first: bool) -> Tree {
        let tips = 2 * cherries + 1;
        let spine = cherries - 1;
        let first_cherry = tips;
        let first_spine = tips + cherries;
        let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); tips + cherries + spine];
        let mut join = |a: usize, b: usize| {
            adj[a].push((b, 0.1));
            adj[b].push((a, 0.1));
        };
        for c in 0..cherries {
            join(first_cherry + c, 1 + 2 * c);
            join(first_cherry + c, 2 + 2 * c);
        }
        join(0, first_spine);
        for s in 0..spine {
            let node = first_spine + s;
            let cherry = first_cherry + s;
            let rest = if s + 1 < spine {
                first_spine + s + 1
            } else {
                first_cherry + cherries - 1
            };
            // the order of the neighbours is the order the children are listed in
            if cherry_first {
                join(node, cherry);
                join(node, rest);
            } else {
                join(node, rest);
                join(node, cherry);
            }
        }
        Tree { adj, tips }
    }

    #[test]
    fn the_larger_subtree_is_done_first_so_that_few_vectors_are_held_whatever_the_order_the_children_are_listed_in(
    ) {
        let q = rate_matrix(2.0, &[0.25; 4]);
        for cherries in [8usize, 100, 400] {
            let bound = (cherries as f64).log2().ceil() as usize + 2;
            for cherry_first in [true, false] {
                let tree = spine_of_cherries(cherries, cherry_first);
                let plan = Plan::new(&tree, &q, &[0.25; 4]).expect("a plan");
                assert_eq!(plan.steps.len(), 2 * cherries - 1);
                let peak = peak_live(&plan);
                assert!(
                    peak <= bound,
                    "{cherries} cherries, cherry first {cherry_first}: {peak} vectors held at once (bound {bound})"
                );
            }
        }
    }

    #[test]
    fn every_node_is_done_after_its_children_and_exactly_once() {
        let q = rate_matrix(2.0, &[0.25; 4]);
        let tree = spine_of_cherries(50, true);
        let plan = Plan::new(&tree, &q, &[0.25; 4]).expect("a plan");
        let mut done = vec![false; plan.nodes];
        for step in &plan.steps {
            assert!(!done[step.node], "node {} twice", step.node);
            for child in step.children.iter().filter(|c| c.tip.is_none()) {
                assert!(done[child.node], "node {} before its child", step.node);
            }
            done[step.node] = true;
        }
        assert!(done[plan.first]);
        assert_eq!(done.iter().filter(|d| **d).count(), plan.steps.len());
    }
}
