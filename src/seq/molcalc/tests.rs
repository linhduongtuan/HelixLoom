//! The calculator against the original. `golden/*.txt` were made by running the original's own procedures (cut out of
//! AppMain.tcl unchanged, run with the system `tclsh` 8.5.9, widgets stood in for) on random inputs; each line holds the
//! inputs and what the procedure left in its variables. The harness is in the notes of MECHANICS.md.

use super::format::{g3, number, tcl_double};
use super::mix::{master_mix, Extra, MixEdit, MixInput};
use super::molarity::{molarity, set_units, Amount, MolarityInput, Per, Unit, Weight};
use super::presets::{DEFAULTS, MIX_PRESETS, MOLAR_PRESETS};
use super::reaction::{molar_reaction, DnaUnits, Fragment, MolarInput};
use super::stock::{update, Changed, Conc, Litres, StockInput};

const MOLAR: &str = include_str!("golden/molar.txt");
const MIX: &str = include_str!("golden/mix.txt");
const MOLARITY: &str = include_str!("golden/molarity.txt");
const STOCK: &str = include_str!("golden/stock.txt");
const UNITS: &str = include_str!("golden/units.txt");
const FORMAT: &str = include_str!("golden/format.txt");

fn fields(line: &str) -> Vec<&str> {
    line.split('\t').collect()
}

#[test]
fn numbers_are_written_as_tcl_writes_them() {
    // each of these was printed by `expr` in tclsh 8.5.9
    for (value, written) in [
        (1e15, "1000000000000000.0"),
        (1e16, "10000000000000000.0"),
        (1e17, "1e+17"),
        (123456789012345678.0, "1.2345678901234568e+17"),
        (0.0001, "0.0001"),
        (0.00001, "1e-5"),
        (0.000123, "0.000123"),
        (100.0, "100.0"),
        (2.5, "2.5"),
        (1e22, "1e+22"),
        (-0.5, "-0.5"),
        (1e-7, "1e-7"),
        (12345.678, "12345.678"),
        (0.1, "0.1"),
        (123_456_789.123_456_79, "123456789.12345679"),
        (1e300, "1e+300"),
        (5e-324, "5e-324"),
        (1.0 / 3.0, "0.3333333333333333"),
        (0.0, "0.0"),
    ] {
        assert_eq!(tcl_double(value), written, "{value:e}");
    }
    assert_eq!(number(" 2.5 "), Some(2.5));
    assert_eq!(number(".5"), Some(0.5));
    assert_eq!(number("5."), Some(5.0));
    for not in ["", " ", ".", "-", "abc", "inf", "nan", "1e", "0x10", "1,5"] {
        assert_eq!(number(not), None, "{not:?}");
    }
}

#[test]
fn percent_g_with_three_digits_is_the_one_of_c() {
    let mut cases = 0;
    for line in FORMAT.lines() {
        let f = fields(line);
        let value = number(f[1]).expect("a number");
        assert_eq!(g3(value), f[2], "{}", f[1]);
        cases += 1;
    }
    assert_eq!(cases, 40);
    for (value, written) in [
        (999.5, "1e+03"),
        (99.95, "100"),
        (0.0, "0"),
        (-1234.5, "-1.23e+03"),
    ] {
        assert_eq!(g3(value), written);
    }
}

#[test]
fn a_prefix_is_chosen_the_way_the_original_chooses_it() {
    let mut cases = 0;
    for line in UNITS.lines() {
        let f = fields(line);
        let value = number(f[1]).expect("a number");
        assert_eq!(
            set_units(value, f[2]).as_deref(),
            Some(f[3]),
            "{} {}",
            f[1],
            f[2]
        );
        cases += 1;
    }
    assert_eq!(cases, 40);
    // the original finds the power of ten with natural logs: 1e-3 comes out a hair above -3, so it is milli
    assert_eq!(set_units(1e-3, "g").as_deref(), Some("1 mg"));
    assert_eq!(set_units(0.0, "mole").as_deref(), Some("0 pmole"));
    assert_eq!(set_units(-1.0, "g"), None);
}

#[test]
fn molar_ratio_reactions_equal_what_the_original_works_out() {
    let mut cases = 0;
    for line in MOLAR.lines() {
        let f = fields(line);
        assert_eq!((f[0], f[8]), ("molar", "ok"), "{line}");
        let fragments: Vec<Fragment> = f[7]
            .split(';')
            .map(|row| {
                let p: Vec<&str> = row.split(',').collect();
                Fragment {
                    name: String::new(),
                    size_kb: p[0].into(),
                    conc: p[1].into(),
                    ratio: p[2].into(),
                }
            })
            .collect();
        let out = molar_reaction(&MolarInput {
            fragments: &fragments,
            reaction_volume: f[1],
            final_amount: f[2],
            units: DnaUnits::ALL
                .into_iter()
                .find(|u| u.label() == f[3])
                .expect("units"),
            min_pipette: f[4],
            enzyme_x: f[5],
            add_water: f[6] == "1",
        });
        for (i, want) in f[9].split(';').enumerate() {
            let w: Vec<&str> = want.split(',').collect();
            let row = &out.rows[i];
            assert_eq!(
                (row.volume.as_str(), row.amount.as_str(), row.ratio.as_str()),
                (w[0], w[1], w[2]),
                "row {i} of {line}"
            );
        }
        let t: Vec<&str> = f[10].split('|').collect();
        assert_eq!(
            (
                out.amount_label,
                out.water.as_str(),
                out.intermediate.as_str(),
                out.dna_volume.as_str(),
                out.enzyme_volume.as_str()
            ),
            (t[0], t[1], t[2], t[3], t[4]),
            "{line}"
        );
        cases += 1;
    }
    assert_eq!(cases, 250);
}

#[test]
fn master_mixes_equal_what_the_original_works_out() {
    let (mut cases, mut skipped) = (0, 0);
    for line in MIX.lines() {
        let f = fields(line);
        assert_eq!((f[0], f[7]), ("mix", "ok"), "{line}");
        let components: Vec<String> = f[6].split(',').map(String::from).collect();
        // the original recalculates the single volume for any trigger that holds "frag_vo"
        let edit = if f[1] == "single_reaction_volume" {
            MixEdit::SingleVolume
        } else if f[1].contains("frag_vo") {
            MixEdit::Component
        } else {
            MixEdit::Other
        };
        let out = master_mix(
            &MixInput {
                components: &components,
                single_volume: f[2],
                reactions: f[3],
                extra: f[4],
                extra_unit: if f[5] == "ul" {
                    Extra::Ul
                } else {
                    Extra::Percent
                },
            },
            edit,
        );
        assert_eq!(out.single_volume, f[8], "{line}");
        // an extra in ul over whole numbers is where the original divides whole numbers and adds nothing: see `mix.rs`
        let whole = |text: &str| !text.is_empty() && text.bytes().all(|b| b.is_ascii_digit());
        if f[5] == "ul" && whole(f[4]) && whole(f[3]) && whole(f[8]) {
            skipped += 1;
            continue;
        }
        for (i, pair) in f[9].split(';').enumerate() {
            let (components, volumes) = pair.split_once(',').expect("a pair");
            assert_eq!(out.components[i], components, "component {i} of {line}");
            assert_eq!(out.volumes[i], volumes, "volume {i} of {line}");
        }
        assert_eq!(out.total, f[10], "{line}");
        cases += 1;
    }
    assert_eq!(cases + skipped, 150);
    assert!(cases > 100, "only {cases} of 150 were compared");
}

#[test]
fn an_extra_in_microlitres_is_added_whatever_way_the_numbers_are_typed() {
    // 12 reactions of 50 ul with 20 ul extra: the original adds nothing for whole numbers; here the 17.99 ul component
    // grows by 20 / 600 of itself
    let components: Vec<String> = ["17.99", "", "", "", "", "", "", "", "", ""]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let out = master_mix(
        &MixInput {
            components: &components,
            single_volume: "50",
            reactions: "12",
            extra: "20",
            extra_unit: Extra::Ul,
        },
        MixEdit::Other,
    );
    assert_eq!(out.volumes[0], "223.08");
    // and as a share of the whole it is the same as the original: 20 % more
    let out = master_mix(
        &MixInput {
            extra_unit: Extra::Percent,
            ..MixInput {
                components: &components,
                single_volume: "50",
                reactions: "12",
                extra: "20",
                extra_unit: Extra::Ul,
            }
        },
        MixEdit::Other,
    );
    assert_eq!(out.volumes[0], "259.06");
}

#[test]
fn conversions_equal_what_the_original_works_out_except_where_it_is_wrong() {
    let (mut cases, mut kept_from_the_original) = (0, 0);
    for line in MOLARITY.lines() {
        let f = fields(line);
        assert_eq!((f[0], f[11]), ("molarity", "ok"), "{line}");
        let find = |label: &str| Unit::FIRST.into_iter().find(|u| u.label() == label);
        let out = molarity(&MolarityInput {
            first: f[1],
            unit1: find(f[2]).expect("a unit"),
            per1: Per::ALL.into_iter().find(|p| p.label() == f[3]),
            unit2: find(f[4]).expect("a unit"),
            per2: Per::ALL
                .into_iter()
                .find(|p| p.label() == f[5])
                .expect("a volume"),
            weight: Weight::ALL
                .into_iter()
                .find(|w| w.label() == f[6])
                .expect("a weight"),
            mol_wt: f[7],
            sequence: f[8],
            amount: f[9],
            amount_unit: Amount::ALL
                .into_iter()
                .find(|a| a.label() == f[10])
                .expect("a unit"),
        });
        assert_eq!(out.second, f[12], "{line}");
        assert_eq!(out.mol_wt, f[13], "{line}");
        // moles per litre alone (key 3) is where the original divides by the volume: see `second_line`
        if out.key != 3 {
            assert_eq!(out.detail, f[14], "{line}");
            kept_from_the_original += 1;
        }
        cases += 1;
    }
    assert_eq!(cases, 300);
    assert!(kept_from_the_original > 250, "{kept_from_the_original}");
}

#[test]
fn moles_per_litre_are_multiplied_by_the_volume_where_the_original_divides() {
    // 1 mmole/L in 0.5 L is 0.5 mmole (500 umole); the original answers "2 mmole"
    let out = molarity(&MolarityInput {
        first: "1",
        unit1: Unit::Mmole,
        per1: Some(Per::L),
        unit2: Unit::Umole,
        per2: Per::L,
        weight: Weight::GramsPerMole,
        mol_wt: "1000",
        sequence: "",
        amount: "0.5",
        amount_unit: Amount::L,
    });
    assert_eq!((out.key, out.second.as_str()), (3, "1e+03"));
    assert_eq!(out.detail, " =  500 umole ");
    // and the litres that hold a number of moles are the moles over the concentration: 1 umole at 1 mmole/L is 1 mL
    let out = molarity(&MolarityInput {
        amount: "1",
        amount_unit: Amount::Umole,
        ..MolarityInput {
            first: "1",
            unit1: Unit::Mmole,
            per1: Some(Per::L),
            unit2: Unit::Umole,
            per2: Per::L,
            weight: Weight::GramsPerMole,
            mol_wt: "1000",
            sequence: "",
            amount: "",
            amount_unit: Amount::L,
        }
    });
    assert_eq!(out.detail, " =  1 mL ");
}

#[test]
fn stock_dilutions_equal_what_the_original_works_out() {
    let conc = |label: &str| {
        Conc::ALL
            .into_iter()
            .find(|c| c.label() == label)
            .expect("a unit")
    };
    let mut cases = 0;
    for line in STOCK.lines() {
        let f = fields(line);
        assert_eq!((f[0], f[8]), ("stock", "ok"), "{line}");
        let changed = match f[1] {
            "stock_units" => Some(Changed::Stock),
            "final_conc_units" => Some(Changed::Final),
            _ => None,
        };
        let (stock_unit, final_unit, out) = update(
            &StockInput {
                stock: f[2],
                stock_unit: conc(f[3]),
                final_conc: f[4],
                final_unit: conc(f[5]),
                final_volume: f[6],
                volume_unit: Litres::ALL
                    .into_iter()
                    .find(|l| l.label() == f[7])
                    .expect("a volume"),
            },
            changed,
        );
        assert_eq!(
            (stock_unit.label(), final_unit.label()),
            (f[9], f[10]),
            "{line}"
        );
        assert_eq!(
            (out.stock.as_str(), out.water.as_str()),
            (f[11], f[12]),
            "{line}"
        );
        cases += 1;
    }
    assert_eq!(cases, 150);
}

#[test]
fn the_reactions_and_first_values_are_the_originals() {
    let names: Vec<&str> = MOLAR_PRESETS.iter().map(|p| p.name).collect();
    assert_eq!(
        names,
        [
            "HiFi 2-3 fragment",
            "HiFi 4-6 fragment",
            "Gibson 2-3 fragment",
            "Gibson 4-6 fragment",
            "SapTrap",
            "Golden Gate",
            "T4 Ligation",
            "NEB Quick Ligation"
        ]
    );
    assert_eq!(DEFAULTS.molar_reaction, MOLAR_PRESETS[0].name);
    assert_eq!(MIX_PRESETS[0].name, DEFAULTS.mix_reaction);
    // Gibson mixes are 1.333x (3:1) and HiFi 2x (1:1); the ligations are in ng/ul
    assert_eq!(MOLAR_PRESETS[2].enzyme_x, "1.333");
    assert_eq!(MOLAR_PRESETS[0].enzyme_x, "2");
    assert!(MOLAR_PRESETS[6..]
        .iter()
        .all(|p| p.units == DnaUnits::NgPerUl));
    // Phusion adds up to the 50 ul that a reaction of it is
    let sum: f64 = MIX_PRESETS[0]
        .components
        .iter()
        .filter_map(|(_, v)| number(v))
        .sum();
    assert_eq!(tcl_double(sum), "50.0");
}

#[test]
fn a_hifi_assembly_of_two_fragments_works_out_by_hand() {
    // 2 fragments of 3 kb and 1 kb at 100 ng/ul, equimolar, 40 fmole in 4 ul with a 2x enzyme mix
    let fragments = [
        Fragment {
            name: "a".into(),
            size_kb: "3".into(),
            conc: "100".into(),
            ratio: "1".into(),
        },
        Fragment {
            name: "b".into(),
            size_kb: "1".into(),
            conc: "100".into(),
            ratio: String::new(),
        },
    ];
    let out = molar_reaction(&MolarInput {
        fragments: &fragments,
        reaction_volume: "4",
        final_amount: "40",
        units: DnaUnits::Fmole,
        min_pipette: "0.3",
        enzyme_x: "2",
        add_water: true,
    });
    // 20 fmole each: 0.02 pmole x (3 x 617.96 + 36.04) ng/pmole is 37.80 ng, 0.378 ul at 100 ng/ul, and the 1 kb one is
    // 0.02 x 654 ng = 0.131 ul. The smaller is under the 0.3 ul that can be pipetted, so everything is scaled by 0.3 /
    // 0.1308 = 2.29: 0.87 and 0.30 ul of 45.87 fmole each, in a DNA mix of 4.59 ul (2.29 x the 2 ul that is half of the
    // reaction) with 3.42 ul of water. These are the numbers the original's own procedure gives for these inputs.
    assert_eq!(out.rows[1].ratio, "1", "a ratio left empty is 1");
    assert_eq!(
        (out.rows[0].volume.as_str(), out.rows[0].amount.as_str()),
        ("0.87", "45.87")
    );
    assert_eq!(
        (out.rows[1].volume.as_str(), out.rows[1].amount.as_str()),
        ("0.30", "45.87")
    );
    assert_eq!(out.amount_label, "fmoles");
    assert_eq!(
        (out.enzyme_volume.as_str(), out.dna_volume.as_str()),
        ("2.0", "2.0")
    );
    assert_eq!(
        (out.water.as_str(), out.intermediate.as_str()),
        ("3.42", "4.59")
    );
}

#[test]
fn the_first_conversion_of_the_window_is_one_milligram_and_one_micromole() {
    // 1 mmole/L = 1 g/L at 1000 g/mol, and 1 mL of it holds 1 mg, which is 1 umole (the original's first screen)
    let out = molarity(&MolarityInput {
        first: DEFAULTS.molarity_first,
        unit1: Unit::Mmole,
        per1: Some(Per::L),
        unit2: Unit::G,
        per2: Per::L,
        weight: Weight::GramsPerMole,
        mol_wt: DEFAULTS.molarity_mol_wt,
        sequence: "",
        amount: DEFAULTS.molarity_amount,
        amount_unit: Amount::Ml,
    });
    assert_eq!((out.second.as_str(), out.key), ("1", 7));
    assert_eq!(out.detail, " =  1 mg , 1 umole ");
}

const CALC: &str = include_str!("golden/calc.txt");

#[test]
fn lines_of_arithmetic_equal_what_the_original_gives() {
    use super::expr::calculate;
    let (mut cases, mut skipped) = (0, 0);
    for line in CALC.lines() {
        let (text, want) = line.split_once('\t').expect("a line");
        // Tcl's whole numbers have no limit; here they stop at 128 bits
        let too_big = want.len() > 38
            && want
                .trim_start_matches('-')
                .bytes()
                .all(|b| b.is_ascii_digit());
        // a number typed as the chosen branch of `?:` in an expression that is not all constants comes back as it was
        // typed (`4.`, `70.E-0`, `2.25e-4`); here it is a number (`4`, `70`, `0.000225`): what is not the way Tcl writes
        // its own number is that
        let typed = number(want).is_some_and(|v| {
            let written = tcl_double(v);
            written.strip_suffix(".0").unwrap_or(&written) != want
        });
        if too_big || typed {
            skipped += 1;
            continue;
        }
        assert_eq!(calculate(text), want, "{text:?}");
        cases += 1;
    }
    assert_eq!(cases + skipped, 3000);
    assert!(cases > 2980, "{cases}");
}

#[test]
fn a_line_too_long_or_nested_too_deeply_is_a_mistake_and_not_a_crash() {
    use super::expr::{calculate, MISTAKE};
    let nest =
        |open: &str, close: &str, n: usize| format!("{}1{}", open.repeat(n), close.repeat(n));
    // what is typed by hand is nested a few levels
    assert_eq!(calculate(&nest("(", ")", 40)), "1");
    assert_eq!(calculate(&nest("sqrt(", ")", 40)), "1");
    // a line of brackets, calls, signs, powers or conditions, as long as a line may be, does not use up the stack
    for line in [
        nest("(", ")", 400),
        nest("sqrt(", ")", 150),
        format!("{}1", "-".repeat(900)),
        format!("{}1", "2**".repeat(300)),
        format!("{}1{}", "1?".repeat(300), ":1".repeat(300)),
    ] {
        assert_eq!(calculate(&line), MISTAKE, "{} characters", line.len());
    }
    // a long sum is not nested; only the length of the line is limited
    assert_eq!(calculate(&format!("{}1", "1+".repeat(400))), "401");
    assert_eq!(calculate(&format!("{}1", "1+".repeat(600))), MISTAKE);
}

#[test]
fn whatever_is_typed_gives_a_result_or_the_mistake_line_and_never_a_crash() {
    use super::expr::calculate;
    // pieces of the kinds of lines a person types, and of the lines that should only give the mistake line
    const PIECES: [&str; 44] = [
        "0",
        "1",
        "7",
        "9",
        "10",
        ".",
        "1.5",
        "2.25e-4",
        "1e",
        "e",
        "E",
        "x",
        "0x",
        "0x1f",
        "$pi",
        "$e",
        "$x",
        " ",
        "(",
        ")",
        ",",
        "+",
        "-",
        "*",
        "/",
        "%",
        "**",
        "<<",
        ">>",
        "&",
        "|",
        "^",
        "~",
        "!",
        "?",
        ":",
        "<=",
        "==",
        "sqrt(",
        "abs(",
        "atan2(",
        "max(",
        "é٣€",
        "170141183460469231731687303715884105727",
    ];
    const BIG: [&str; 4] = [
        "170141183460469231731687303715884105728",
        "0x7fffffffffffffffffffffffffffffff",
        "-9223372036854775808",
        "99999999999999999999",
    ];
    // a small generator of its own, so that the lines are the same every time
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let mut next = move |n: usize| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % n as u64) as usize
    };
    for _ in 0..30_000 {
        let line: String = (0..1 + next(14))
            .map(|_| match next(12) {
                0 => BIG[next(BIG.len())],
                _ => PIECES[next(PIECES.len())],
            })
            .collect();
        assert!(!calculate(&line).is_empty(), "{line:?}");
    }
}

/// A small generator of its own, so that the same fields are tried every time.
struct Dice(u64);

impl Dice {
    fn below(&mut self, n: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % n as u64) as usize
    }

    fn pick<T: Copy>(&mut self, items: &[T]) -> T {
        items[self.below(items.len())]
    }
}

#[test]
fn whatever_is_typed_in_the_fields_gives_numbers_or_nothing_and_never_a_crash() {
    use super::mix::{master_mix, Extra, MixEdit, MixInput, COMPONENTS};
    use super::molarity::{molarity, Amount, MolarityInput, Per, Unit, Weight};
    use super::reaction::{molar_reaction, DnaUnits, Fragment, MolarInput, FRAGMENTS};
    use super::stock::{stock_dilution, update, Changed, Conc, Litres, StockInput};
    const TEXT: [&str; 20] = [
        "",
        " ",
        "0",
        "-1",
        "1",
        "2.5",
        "10",
        "1e999",
        "-1e999",
        "1e-999",
        "abc",
        ".",
        "-",
        "1e",
        "NaN",
        "inf",
        "0.0000001",
        "1000000",
        "3",
        "100",
    ];
    const SEQUENCES: [&str; 5] = ["", "ACGT", "acgtu", "ACGTXYZ 123", "UUUUTTTT"];
    let mut dice = Dice(0x9e37_79b9_7f4a_7c15);
    for _ in 0..20_000 {
        let fragments: Vec<Fragment> = (0..1 + dice.below(FRAGMENTS + 2))
            .map(|_| Fragment {
                name: "f".into(),
                size_kb: dice.pick(&TEXT).into(),
                conc: dice.pick(&TEXT).into(),
                ratio: dice.pick(&TEXT).into(),
            })
            .collect();
        let reaction = molar_reaction(&MolarInput {
            fragments: &fragments,
            reaction_volume: dice.pick(&TEXT),
            final_amount: dice.pick(&TEXT),
            units: dice.pick(&DnaUnits::ALL),
            min_pipette: dice.pick(&TEXT),
            enzyme_x: dice.pick(&TEXT),
            add_water: dice.below(2) == 0,
        });
        assert_eq!(reaction.rows.len(), fragments.len().min(FRAGMENTS));

        let components: Vec<String> = (0..1 + dice.below(COMPONENTS + 2))
            .map(|_| dice.pick(&TEXT).into())
            .collect();
        let edit = dice.pick(&[MixEdit::Component, MixEdit::SingleVolume, MixEdit::Other]);
        let mix = master_mix(
            &MixInput {
                components: &components,
                single_volume: dice.pick(&TEXT),
                reactions: dice.pick(&TEXT),
                extra: dice.pick(&TEXT),
                extra_unit: dice.pick(&[Extra::Percent, Extra::Ul]),
            },
            edit,
        );
        assert_eq!(mix.components.len(), components.len().min(COMPONENTS));

        let per1 = dice.pick(&[None, Some(Per::L), Some(Per::Ml), Some(Per::Ul)]);
        molarity(&MolarityInput {
            first: dice.pick(&TEXT),
            unit1: dice.pick(&Unit::FIRST),
            per1,
            unit2: dice.pick(&Unit::SECOND),
            per2: dice.pick(&Per::ALL),
            weight: dice.pick(&Weight::ALL),
            mol_wt: dice.pick(&TEXT),
            sequence: dice.pick(&SEQUENCES),
            amount: dice.pick(&TEXT),
            amount_unit: dice.pick(&Amount::ALL),
        });

        let dilution = StockInput {
            stock: dice.pick(&TEXT),
            stock_unit: dice.pick(&Conc::ALL),
            final_conc: dice.pick(&TEXT),
            final_unit: dice.pick(&Conc::ALL),
            final_volume: dice.pick(&TEXT),
            volume_unit: dice.pick(&Litres::ALL),
        };
        stock_dilution(&dilution);
        update(
            &dilution,
            dice.pick(&[None, Some(Changed::Stock), Some(Changed::Final)]),
        );
    }
}

#[test]
fn zero_to_a_negative_power_is_a_mistake_with_the_operator_and_infinity_with_the_function() {
    use super::expr::{calculate, MISTAKE};
    // checked with the original's own `calc_update` under tclsh: "exponentiation of zero by negative power"
    for line in ["0**-1", "0 ** -2", "-0**-1", "0.0**-1", "0**(0-3)"] {
        assert_eq!(calculate(line), MISTAKE, "{line}");
    }
    assert_eq!(calculate("pow(0,-1)"), "Inf");
    assert_eq!(calculate("0.5**-1"), "2");
    assert_eq!(calculate("0**0"), "1");
}
