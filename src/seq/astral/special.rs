//! The numbers a branch gets from its quartet frequencies: a length in coalescent units and the local posterior probability
//! of each of the three ways to resolve the branch (Sayyari and Mirarab 2016, Molecular Biology and Evolution 33: 1654).
//!
//! The gene trees show the three quartet topologies in the frequencies `f1 f2 f3` (of `en = f1 + f2 + f3` genes). Under the
//! coalescent a gene tree has the branch's own topology with probability `p`, and each of the others with `(1 - p) / 2`,
//! where `p = 1 - 2/3 e^-t` lies between 1/3 (no branch) and 1 (a long one). With every topology equally likely beforehand
//! and `p` uniform on 1/3 to 1, the likelihood of the frequencies for the topology `i` is
//! `3/2 * 2^-(en - f_i) * B(f_i + 1, en - f_i + 1) * (1 - I_{1/3}(f_i + 1, en - f_i + 1))`, with `B` the beta function and
//! `I` the regularized incomplete beta function; the posteriors are these three, scaled to add up to 1.

use std::f64::consts::PI;

/// The width of the Lanczos approximation of the log-gamma function used below.
const LANCZOS_G: f64 = 7.0;
const LANCZOS: [f64; 9] = [
    0.999_999_999_999_809_9,
    676.520_368_121_885_1,
    -1_259.139_216_722_402_8,
    771.323_428_777_653_1,
    -176.615_029_162_140_6,
    12.507_343_278_686_905,
    -0.138_571_095_265_720_12,
    9.984_369_578_019_572e-6,
    1.505_632_735_149_311_6e-7,
];
/// How close two steps of the continued fraction must be, and how many it may take.
const FRACTION_EPSILON: f64 = 1e-16;
const FRACTION_STEPS: usize = 20_000;
const TINY: f64 = 1e-300;

/// The natural logarithm of the gamma function, for positive `x`.
pub fn ln_gamma(x: f64) -> f64 {
    if x < 0.5 {
        return PI.ln() - (PI * x).sin().ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let t = x + LANCZOS_G + 0.5;
    let series = LANCZOS[0]
        + LANCZOS[1..]
            .iter()
            .enumerate()
            .map(|(i, c)| c / (x + (i + 1) as f64))
            .sum::<f64>();
    0.5 * (2.0 * PI).ln() + (x + 0.5) * t.ln() - t + series.ln()
}

fn ln_beta(a: f64, b: f64) -> f64 {
    ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b)
}

/// The continued fraction of the incomplete beta function, by Lentz's method.
fn beta_fraction(a: f64, b: f64, x: f64) -> f64 {
    let (sum, plus, minus) = (a + b, a + 1.0, a - 1.0);
    let guard = |v: f64| if v.abs() < TINY { TINY } else { v };
    let mut c = 1.0;
    let mut d = 1.0 / guard(1.0 - sum * x / plus);
    let mut h = d;
    for m in 1..=FRACTION_STEPS {
        let m = m as f64;
        let twice = 2.0 * m;
        let even = m * (b - m) * x / ((minus + twice) * (a + twice));
        d = 1.0 / guard(1.0 + even * d);
        c = guard(1.0 + even / c);
        h *= d * c;
        let odd = -(a + m) * (sum + m) * x / ((a + twice) * (plus + twice));
        d = 1.0 / guard(1.0 + odd * d);
        c = guard(1.0 + odd / c);
        let step = d * c;
        h *= step;
        if (step - 1.0).abs() < FRACTION_EPSILON {
            break;
        }
    }
    h
}

/// The natural logarithm of the regularized incomplete beta function `I_x(a, b)`, for `0 < x < 1`.
pub fn ln_beta_i(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if x >= 1.0 {
        return 0.0;
    }
    let front = a * x.ln() + b * (1.0 - x).ln() - ln_beta(a, b);
    if x < (a + 1.0) / (a + b + 2.0) {
        front + (beta_fraction(a, b, x) / a).ln()
    } else {
        // I_x(a, b) = 1 - I_{1-x}(b, a): the fraction converges faster on the other side
        let other = (front + (beta_fraction(b, a, 1.0 - x) / b).ln()).exp();
        (-other).ln_1p()
    }
}

/// The local posterior probabilities of the three topologies of a branch from its quartet frequencies.
pub fn posteriors(f: [f64; 3]) -> [f64; 3] {
    let en: f64 = f.iter().sum();
    let logs = f.map(|fi| {
        let (a, b) = (fi + 1.0, en - fi + 1.0);
        // ln(1 - I_{1/3}(a, b)) is ln I_{2/3}(b, a)
        1.5_f64.ln() - (en - fi) * std::f64::consts::LN_2
            + ln_beta(a, b)
            + ln_beta_i(b, a, 2.0 / 3.0)
    });
    let most = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let weights = logs.map(|l| (l - most).exp());
    let total: f64 = weights.iter().sum();
    weights.map(|w| w / total)
}

/// The length of a branch in coalescent units from the frequency `f1` of its own topology among `en` genes: the value that
/// is most probable under a Yule prior of rate `lambda` (ASTRAL's default is 0.5; 0 gives the maximum likelihood length,
/// which is infinite when every gene agrees). A branch whose topology is not the most frequent one has no length.
pub fn coalescent_length(f1: f64, en: f64, lambda: f64) -> f64 {
    let agree = f1 / (en + 2.0 * lambda);
    let shortfall = 1.5 * (1.0 - agree);
    if shortfall >= 1.0 {
        0.0
    } else if shortfall <= 0.0 {
        f64::INFINITY
    } else {
        -shortfall.ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, relative: f64) -> bool {
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(f64::MIN_POSITIVE)
    }

    #[test]
    fn the_log_gamma_function_gives_factorials_and_the_square_root_of_pi() {
        assert!(close(ln_gamma(5.0), 24.0_f64.ln(), 1e-13));
        assert!(close(ln_gamma(0.5), PI.sqrt().ln(), 1e-13));
        assert!(close(ln_gamma(101.0), 363.739_375_555_563_5, 1e-13));
        assert!(close(ln_gamma(0.1), 2.252_712_651_734_206, 1e-13));
    }

    #[test]
    fn the_incomplete_beta_function_matches_known_values_in_both_tails() {
        // I_x(1, 1) = x; I_x(2, 1) = x^2; I_x(1, 2) = 1 - (1 - x)^2
        assert!(close(ln_beta_i(1.0, 1.0, 0.3).exp(), 0.3, 1e-13));
        assert!(close(ln_beta_i(2.0, 1.0, 0.3).exp(), 0.09, 1e-13));
        assert!(close(ln_beta_i(1.0, 2.0, 0.7).exp(), 0.91, 1e-13));
        // the symmetric case: I_{1/2}(a, a) = 1/2
        assert!(close(ln_beta_i(40.0, 40.0, 0.5).exp(), 0.5, 1e-12));
        // far in the lower tail: a hand-checkable bound, I_x(a, b) < x^a / (a B(a, b)) * (1 - x)^b / (1 - x)
        let tail = ln_beta_i(300.0, 100.0, 0.3);
        assert!(tail < -40.0 && tail.is_finite(), "{tail}");
    }

    #[test]
    fn posteriors_equal_the_ones_astral_prints() {
        // one gene that agrees: 2/3 for the branch, 1/6 for each of the others
        let one = posteriors([1.0, 0.0, 0.0]);
        assert!(close(one[0], 2.0 / 3.0, 1e-13) && close(one[1], 1.0 / 6.0, 1e-13));
        // frequencies of an inner branch of the primates data and what ASTRAL 5.7.8 printed for them
        let pp = posteriors([164.861_111_111_111_14, 150.75, 108.388_888_888_888_89]);
        assert!(close(pp[0], 0.9247407068234015, 1e-9), "{pp:?}");
        assert!(close(pp[1], 0.06959483623204127, 1e-9), "{pp:?}");
        assert!(close(pp[2], 0.00566445694455617, 1e-9), "{pp:?}");
        // a branch that is nearly certain: the others are tiny. ASTRAL prints 2.28e-54 and 0 for them, which is what
        // subtracting from 1 in doubles leaves; the exact values are about 2.5e-54 and 2.3e-54
        let sure = posteriors([298.0, 67.090_909_090_909_07, 58.909_090_909_090_94]);
        assert!(
            sure[0] > 1.0 - 1e-12 && sure[1] < 1e-50 && sure[2] < 1e-50,
            "{sure:?}"
        );
        assert!(sure[1] > 1e-56 && sure[2] > 1e-56, "{sure:?}");
        assert!((sure.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        // no gene says anything: a third each
        assert!(posteriors([0.0; 3])
            .iter()
            .all(|&p| close(p, 1.0 / 3.0, 1e-13)));
    }

    #[test]
    fn posteriors_equal_the_formula_computed_to_forty_digits() {
        // `(frequencies, posteriors)` by the formula of Sayyari and Mirarab (2016) in 40-digit arithmetic (Python's mpmath).
        // The first is a branch where ASTRAL 5.7.8 printed 1.4133151e-5 for the smallest: that is 8e-8 off the formula,
        // and this program's value is within 1e-12 of it
        let cases = [
            "7.991285392367723 40.53577938678393 54.84315297685183 | 1.4133150169044268e-5 0.000404849969773153 0.9995810168800578",
            "14 10 11 | 0.599160687031838 0.17660613431833608 0.22423317864982595",
            "35.5 0.75 0.75 | 0.9999999999999986 7.436317710504098e-16 7.436317710504098e-16",
            "22.291666666666664 4.916666666666666 7.791666666666666 | 0.9995449622870919 0.00018201560140197192 0.0002730221115061245",
            "0 3 1 | 0.057971014492753624 0.8405797101449275 0.10144927536231885",
            "1 0 0 | 0.6666666666666666 0.16666666666666666 0.16666666666666666",
            "150.5 149.5 150 | 0.34679277579798906 0.3201509271537959 0.33305629704821504",
            "2 90 3 | 6.2078652967241805e-37 1 6.416540114343792e-37",
        ];
        let numbers = |text: &str| -> Vec<f64> {
            text.split_whitespace()
                .map(|x| x.parse().expect("a number"))
                .collect()
        };
        for case in cases {
            let (f, exact) = case.split_once(" | ").expect("two halves");
            let (f, exact) = (numbers(f), numbers(exact));
            let mine = posteriors([f[0], f[1], f[2]]);
            for (m, e) in mine.iter().zip(&exact) {
                assert!(close(*m, *e, 1e-10), "{f:?}: {mine:?} against {exact:?}");
            }
        }
    }

    #[test]
    fn lengths_follow_the_yule_prior_and_stop_at_zero() {
        // numbers ASTRAL printed for the branches of the primates tree (lambda 0.5)
        assert!(close(
            coalescent_length(399.342_857_142_857_16, 424.0, 0.5),
            2.401_802_054_003_467_5,
            1e-12
        ));
        assert!(close(
            coalescent_length(298.0, 424.0, 0.5),
            0.802_436_974_357_661_1,
            1e-12
        ));
        assert!(close(
            coalescent_length(1.0, 1.0, 0.5),
            0.287_682_072_451_780_9,
            1e-12
        ));
        // a topology that is no more frequent than chance has no length; with no prior, a branch every gene shows is infinite
        assert_eq!(coalescent_length(100.0, 424.0, 0.5), 0.0);
        assert!(coalescent_length(424.0, 424.0, 0.0).is_infinite());
    }
}
