//! A set of taxa as a row of bits: the clusters of the species-tree search and the branches of the gene trees.

/// A set of the taxa `0..n`, kept in 64-bit words.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Bits(Box<[u64]>);

impl Bits {
    /// The empty set of `n` taxa.
    pub fn new(n: usize) -> Bits {
        Bits(vec![0; n.div_ceil(64).max(1)].into_boxed_slice())
    }

    pub fn of(n: usize, items: impl IntoIterator<Item = usize>) -> Bits {
        let mut set = Bits::new(n);
        for i in items {
            set.insert(i);
        }
        set
    }

    pub fn insert(&mut self, i: usize) {
        self.0[i / 64] |= 1 << (i % 64);
    }

    pub fn contains(&self, i: usize) -> bool {
        self.0.get(i / 64).is_some_and(|w| w >> (i % 64) & 1 == 1)
    }

    pub fn count(&self) -> u32 {
        self.0.iter().map(|w| w.count_ones()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.0.iter().all(|&w| w == 0)
    }

    /// How many taxa the two sets share.
    pub fn and_count(&self, other: &Bits) -> u32 {
        self.0
            .iter()
            .zip(other.0.iter())
            .map(|(a, b)| (a & b).count_ones())
            .sum()
    }

    pub fn union(&self, other: &Bits) -> Bits {
        Bits(
            self.0
                .iter()
                .zip(other.0.iter())
                .map(|(a, b)| a | b)
                .collect(),
        )
    }

    /// The taxa of `self` that `other` does not have.
    pub fn minus(&self, other: &Bits) -> Bits {
        Bits(
            self.0
                .iter()
                .zip(other.0.iter())
                .map(|(a, b)| a & !b)
                .collect(),
        )
    }

    pub fn is_disjoint(&self, other: &Bits) -> bool {
        self.0.iter().zip(other.0.iter()).all(|(a, b)| a & b == 0)
    }

    pub fn is_subset_of(&self, other: &Bits) -> bool {
        self.0.iter().zip(other.0.iter()).all(|(a, b)| a & !b == 0)
    }

    /// Every taxon of `0..n` that is not in the set.
    pub fn complement(&self, n: usize) -> Bits {
        Bits::of(n, (0..n).filter(|&i| !self.contains(i)))
    }

    pub fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.0.len() * 64).filter(|&i| self.contains(i))
    }

    pub fn first(&self) -> Option<usize> {
        self.iter().next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sets_count_combine_and_compare_across_word_boundaries() {
        let a = Bits::of(130, [0, 63, 64, 129]);
        let b = Bits::of(130, [63, 64, 100]);
        assert_eq!((a.count(), b.count()), (4, 3));
        assert_eq!(a.and_count(&b), 2);
        assert_eq!(a.union(&b).count(), 5);
        assert_eq!(a.minus(&b).iter().collect::<Vec<_>>(), vec![0, 129]);
        assert!(!a.is_disjoint(&b));
        assert!(Bits::of(130, [63]).is_subset_of(&a));
        assert!(!b.is_subset_of(&a));
        assert_eq!(a.complement(130).count(), 126);
        assert_eq!(a.first(), Some(0));
        assert!(Bits::new(5).is_empty() && Bits::new(0).is_empty());
        assert!(a.contains(129) && !a.contains(130) && !a.contains(5000));
    }
}
