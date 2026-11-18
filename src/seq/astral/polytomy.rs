//! The test for a polytomy (Sayyari and Mirarab 2018, *Genes* 9: 132): whether an inner branch of the species tree is so short
//! that it should be a polytomy.
//!
//! With no branch at all the gene trees show the three ways to resolve the four groups around it equally often. So the three
//! counts of the branch (`weight::resolved_around`: every gene tree gives one in all, divided among the pairings by its
//! resolved quartets) are tested against a third each with a chi-squared statistic, `sum (n_i - n/3)^2 / (n/3)` over the
//! three, on 2 degrees of freedom, where `n` is the number of gene trees that count. The paper uses the test only from 10 gene
//! trees on. A small p-value says the branch is real; a large one does not say it is not (the test may only be too weak).
//! The branches are tested one by one, with no correction for the number of them.

use super::bits::Bits;
use super::genes::Gene;
use super::weight::resolved_around;

/// The fewest gene trees with a resolved quartet around a branch for which it is tested.
pub const MIN_GENES: usize = 10;
/// The level of significance the paper uses: a branch with a p-value under it is taken to have a length.
pub const ALPHA: f64 = 0.05;

/// The outcome of the test on one branch.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Polytomy {
    /// The gene trees that count.
    pub genes: usize,
    pub chi2: f64,
    pub p: f64,
}

impl Polytomy {
    /// Whether the test rejects a polytomy at the level `ALPHA`.
    pub fn rejected(&self) -> bool {
        self.p < ALPHA
    }
}

/// The test on the three counts of `genes` gene trees; none when there are fewer than `MIN_GENES`.
pub fn chi_square(counts: [f64; 3], genes: usize) -> Option<Polytomy> {
    if genes < MIN_GENES {
        return None;
    }
    let expected = counts.iter().sum::<f64>() / 3.0;
    let chi2: f64 = counts.iter().map(|c| (c - expected).powi(2)).sum::<f64>() / expected;
    // the upper tail of the chi-squared distribution on 2 degrees of freedom is exactly exp(-x/2)
    Some(Polytomy {
        genes,
        chi2,
        p: (-chi2 / 2.0).exp(),
    })
}

/// The test for the branch with the four groups `l r | s o` around it.
pub fn around(genes: &[Gene], groups: [&Bits; 4]) -> Option<Polytomy> {
    let [l, r, s, o] = groups;
    let found = resolved_around(genes, l, r, s, o);
    chi_square(found.n, found.genes)
}

/// A p-value as a branch label: three decimals, or the exponent form when that would show nothing (`3.2E-17`).
pub fn p_text(p: f64) -> String {
    if p >= 0.001 {
        let text = format!("{p:.3}");
        text.trim_end_matches('0').trim_end_matches('.').to_string()
    } else if p == 0.0 {
        "0".into()
    } else {
        format!("{p:.1e}").replace('e', "E")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 * a.abs().max(b.abs()).max(1e-300)
    }

    #[test]
    fn equal_counts_cannot_be_told_from_a_polytomy() {
        let t = chi_square([10.0, 10.0, 10.0], 30).expect("tested");
        assert_eq!((t.chi2, t.p, t.genes), (0.0, 1.0, 30));
        assert!(!t.rejected());
    }

    #[test]
    fn the_statistic_and_its_tail_are_the_ones_of_a_chi_squared_on_two_degrees() {
        // 20 5 5 of 30: expected 10 each, (100 + 25 + 25) / 10 = 15, and exp(-7.5) = 5.530843701478336e-4
        let t = chi_square([20.0, 5.0, 5.0], 30).expect("tested");
        assert!(close(t.chi2, 15.0), "{}", t.chi2);
        assert!(close(t.p, 5.530_843_701_478_336e-4), "{}", t.p);
        assert!(t.rejected());
        // the 5% point of a chi-squared on 2 degrees is 5.991464547107979 (a table): counts 10 + x, 10 - x and 10 make it
        // with x^2 / 5, so the p-value there is 0.05
        let x = (5.0 * 5.991_464_547_107_979_f64).sqrt();
        let t = chi_square([10.0 + x, 10.0 - x, 10.0], 30).expect("tested");
        assert!((t.p - 0.05).abs() < 1e-12, "{}", t.p);
        // 3 sides of 12 12 6 (n = 30, expected 10): chi2 = (4 + 4 + 16) / 10 = 2.4, p = exp(-1.2) = 0.30119421191220214
        let t = chi_square([12.0, 12.0, 6.0], 30).expect("tested");
        assert!(close(t.chi2, 2.4) && close(t.p, 0.301_194_211_912_202_14));
        assert!(!t.rejected());
        // the order of the counts does not matter
        assert_eq!(
            chi_square([6.0, 12.0, 12.0], 30).map(|t| t.chi2),
            Some(t.chi2)
        );
    }

    #[test]
    fn a_branch_is_rejected_as_a_polytomy_only_with_a_p_value_under_the_level() {
        let at = |p| Polytomy {
            genes: 30,
            chi2: 0.0,
            p,
        };
        assert!(at(0.049_999).rejected());
        assert!(!at(ALPHA).rejected(), "exactly 0.05 is not under 0.05");
        assert!(!at(0.050_001).rejected());
    }

    #[test]
    fn fewer_than_ten_gene_trees_are_not_tested() {
        assert_eq!(chi_square([9.0, 0.0, 0.0], 9), None);
        assert_eq!(chi_square([0.0, 0.0, 0.0], 0), None);
        let at_ten = chi_square([10.0, 0.0, 0.0], 10).expect("ten are enough");
        // 10 0 0: expected 10/3, (400/9 + 100/9 + 100/9) / (10/3) = 20, p = exp(-10)
        assert!(close(at_ten.chi2, 20.0) && close(at_ten.p, (-10.0_f64).exp()));
    }

    #[test]
    fn a_p_value_is_written_with_three_decimals_or_as_a_power() {
        assert_eq!(p_text(1.0), "1");
        assert_eq!(p_text(0.05), "0.05");
        assert_eq!(p_text(0.30119), "0.301");
        assert_eq!(p_text(0.001), "0.001");
        assert_eq!(p_text(5.53e-4), "5.5E-4");
        assert_eq!(p_text(3.2e-17), "3.2E-17");
        assert_eq!(p_text(0.0), "0");
    }
}
