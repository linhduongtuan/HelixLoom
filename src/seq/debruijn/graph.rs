//! The graph of segments. A segment is an unbranched run of k-mers written out as bases; the graph is *bidirected*: a segment is
//! met in either direction, and a link from the end of one to the start of the next (their last and first k-1 bases are the same) is
//! also a link between the reverse complements, the other way round. A **handle** is a segment in one direction: twice the segment's
//! number, plus one for the reverse.

use crate::seq::dna::revcom;

pub type Handle = usize;

pub fn rc(h: Handle) -> Handle {
    h ^ 1
}

pub fn seg(h: Handle) -> usize {
    h >> 1
}

pub fn handle(id: usize, reverse: bool) -> Handle {
    id * 2 + usize::from(reverse)
}

#[derive(Clone, Debug, PartialEq)]
pub struct Segment {
    /// The bases, at least k of them.
    pub seq: Vec<u8>,
    /// The sum of the counts of the k-mers of the segment: its average is the segment's coverage.
    pub weight: f64,
}

#[derive(Clone, Debug)]
pub struct Graph {
    k: usize,
    segs: Vec<Segment>,
    alive: Vec<bool>,
    /// The handles that follow each handle.
    succ: Vec<Vec<Handle>>,
}

impl Graph {
    pub fn new(k: usize) -> Graph {
        Graph {
            k,
            segs: Vec::new(),
            alive: Vec::new(),
            succ: Vec::new(),
        }
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub fn add(&mut self, seq: Vec<u8>, weight: f64) -> usize {
        self.segs.push(Segment { seq, weight });
        self.alive.push(true);
        self.succ.push(Vec::new());
        self.succ.push(Vec::new());
        self.segs.len() - 1
    }

    /// The numbers of the segments that are there, in order.
    pub fn ids(&self) -> Vec<usize> {
        (0..self.segs.len()).filter(|&i| self.alive[i]).collect()
    }

    /// How many segments are there.
    pub fn count(&self) -> usize {
        self.alive.iter().filter(|&&a| a).count()
    }

    /// How many handles the numbering has room for (some may be of removed segments).
    pub fn handles(&self) -> usize {
        self.segs.len() * 2
    }

    pub fn alive(&self, id: usize) -> bool {
        self.alive.get(id).copied().unwrap_or(false)
    }

    pub fn segment(&self, id: usize) -> &Segment {
        &self.segs[id]
    }

    pub fn length(&self, id: usize) -> usize {
        self.segs[id].seq.len()
    }

    /// The k-mers of the segment.
    pub fn kmers(&self, id: usize) -> usize {
        self.length(id) + 1 - self.k
    }

    /// The average number of times a k-mer of the segment was seen.
    pub fn coverage(&self, id: usize) -> f64 {
        self.segs[id].weight / self.kmers(id) as f64
    }

    pub fn add_weight(&mut self, id: usize, weight: f64) {
        self.segs[id].weight += weight;
    }

    pub fn total_length(&self) -> usize {
        self.ids().into_iter().map(|i| self.length(i)).sum()
    }

    /// The sum of the weights of the segments that are there.
    #[cfg(test)]
    pub fn weight_total(&self) -> f64 {
        self.ids().into_iter().map(|i| self.segs[i].weight).sum()
    }

    /// The bases as they read going the way the handle goes.
    pub fn oriented(&self, h: Handle) -> Vec<u8> {
        let seq = &self.segs[seg(h)].seq;
        if h & 1 == 0 {
            seq.clone()
        } else {
            revcom(seq)
        }
    }

    pub fn succ(&self, h: Handle) -> &[Handle] {
        &self.succ[h]
    }

    /// The handles that lead into `h`.
    pub fn pred(&self, h: Handle) -> Vec<Handle> {
        self.succ[rc(h)].iter().map(|&x| rc(x)).collect()
    }

    pub fn outdeg(&self, h: Handle) -> usize {
        self.succ[h].len()
    }

    pub fn indeg(&self, h: Handle) -> usize {
        self.succ[rc(h)].len()
    }

    /// A link from `a` to `b`, and with it the link from the reverse complement of `b` to that of `a`; a link that is there is not added twice.
    pub fn add_link(&mut self, a: Handle, b: Handle) {
        if !self.succ[a].contains(&b) {
            self.succ[a].push(b);
        }
        if !self.succ[rc(b)].contains(&rc(a)) {
            self.succ[rc(b)].push(rc(a));
        }
    }

    /// Take one link out (and the link between the reverse complements that goes with it); the segments stay.
    pub fn remove_link(&mut self, a: Handle, b: Handle) {
        self.succ[a].retain(|&x| x != b);
        self.succ[rc(b)].retain(|&x| x != rc(a));
    }

    /// Take a segment out, with every link to and from it.
    pub fn remove(&mut self, id: usize) {
        for h in [handle(id, false), handle(id, true)] {
            for x in std::mem::take(&mut self.succ[h]) {
                self.succ[rc(x)].retain(|&y| y != rc(h));
            }
            for x in std::mem::take(&mut self.succ[rc(h)]) {
                self.succ[rc(x)].retain(|&y| y != h);
            }
        }
        self.alive[id] = false;
        self.segs[id].seq = Vec::new();
    }

    /// Every link once: `a` to `b` and the reverse complement of `b` to that of `a` are one link, written as the smaller of the two.
    pub fn links(&self) -> Vec<(Handle, Handle)> {
        let mut out = Vec::new();
        for id in self.ids() {
            for a in [handle(id, false), handle(id, true)] {
                for &b in &self.succ[a] {
                    if (a, b) <= (rc(b), rc(a)) {
                        out.push((a, b));
                    }
                }
            }
        }
        out
    }

    /// A segment that is a whole circle: its end leads to its own start, and to nothing else.
    pub fn circular(&self, id: usize) -> bool {
        let (f, r) = (handle(id, false), handle(id, true));
        self.succ[f] == [f] && self.succ[r] == [r]
    }

    /// Where `a` leads to one segment only, that segment is reached from `a` only, and `joinable` (asked of the two segments' numbers) says
    /// they may be joined: the two are one run.
    fn next(&self, a: Handle, joinable: &dyn Fn(usize, usize) -> bool) -> Option<Handle> {
        let b = *self.succ[a].first().filter(|_| self.succ[a].len() == 1)?;
        (self.succ[rc(b)].len() == 1 && seg(a) != seg(b) && joinable(seg(a), seg(b))).then_some(b)
    }

    /// The graph with each run of segments that has no branch in it made one segment (the numbers change). The bases of a run are those of
    /// its segments, each joined to the one before without the k-1 that they share; the weights add.
    pub fn compact(&self) -> Graph {
        self.compact_where(&|_, _| true)
    }

    /// `compact`, joining two segments only where `joinable` allows it (asked of their numbers, so it must not depend on the direction).
    pub fn compact_where(&self, joinable: &dyn Fn(usize, usize) -> bool) -> Graph {
        let n = self.segs.len();
        let mut place: Vec<Option<(usize, usize)>> = vec![None; n];
        let mut chains: Vec<Vec<Handle>> = Vec::new();
        for s in (0..n).filter(|&s| self.alive[s]) {
            if place[s].is_some() {
                continue;
            }
            // back to where the run starts (a run that is a ring has no start: any place will do)
            let mut start = handle(s, false);
            while let Some(p) = self.next(rc(start), joinable).map(rc) {
                if seg(p) == s {
                    break;
                }
                start = p;
            }
            let chain_no = chains.len();
            let mut chain = vec![start];
            place[seg(start)] = Some((chain_no, 0));
            let mut cur = start;
            while let Some(b) = self.next(cur, joinable) {
                if place[seg(b)].is_some() {
                    break;
                }
                place[seg(b)] = Some((chain_no, chain.len()));
                chain.push(b);
                cur = b;
            }
            chains.push(chain);
        }
        let mut out = Graph::new(self.k);
        for chain in &chains {
            let mut seq = self.oriented(chain[0]);
            for &h in &chain[1..] {
                seq.extend_from_slice(&self.oriented(h)[self.k - 1..]);
            }
            let weight = chain.iter().map(|&h| self.segs[seg(h)].weight).sum();
            out.add(seq, weight);
        }
        let new_handle = |x: Handle| -> Handle {
            let (c, i) = place[seg(x)].unwrap_or((0, 0));
            handle(c, chains[c][i] != x)
        };
        for id in self.ids() {
            for a in [handle(id, false), handle(id, true)] {
                for &b in &self.succ[a] {
                    let ((ca, ia), (cb, ib)) = (
                        place[seg(a)].unwrap_or((0, 0)),
                        place[seg(b)].unwrap_or((0, 0)),
                    );
                    let inside = ca == cb
                        && ((chains[ca][ia] == a && chains[cb][ib] == b && ib == ia + 1)
                            || (chains[ca][ia] == rc(a)
                                && chains[cb][ib] == rc(b)
                                && ia == ib + 1));
                    if !inside {
                        out.add_link(new_handle(a), new_handle(b));
                    }
                }
            }
        }
        out
    }
}

#[cfg(test)]
mod tests;
