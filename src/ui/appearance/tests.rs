use super::*;
use crate::ui::design::{base_fill, residue_fill, trace_color, READABLE};

fn is_default(a: Appearance) -> bool {
    a == Appearance::default()
}

// ---- the settings and the file ---------------------------------------------------------------------------------------------------

#[test]
fn the_defaults_change_nothing_about_the_look_of_either_theme() {
    for dark in [false, true] {
        assert_eq!(Appearance::default().look(dark), design::look(dark));
        let palette = crate::ui::theme::palette(dark);
        let tuned = Appearance::default().palette(palette);
        assert_eq!(tuned.dim, palette.dim);
        assert_eq!(tuned.line, palette.line);
        assert_eq!(tuned.fg, palette.fg);
    }
    assert_eq!(Appearance::default().zoom(), 1.0);
    assert_eq!(
        Appearance::default().density.spacing().interact_height,
        26.0
    );
}

#[test]
fn every_setting_is_written_and_read_back() {
    for theme in ThemeMode::ALL {
        for &scale in &SCALES {
            for &density in Density::ALL {
                for &accent in Accent::ALL {
                    for &contrast in Contrast::ALL {
                        for &colours in Colours::ALL {
                            let a = Appearance {
                                scale,
                                density,
                                accent,
                                contrast,
                                colours,
                            };
                            assert_eq!(Appearance::decode(&a.encode(theme)), (theme, a));
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn a_file_is_written_a_setting_to_a_line_in_words() {
    let a = Appearance {
        scale: 125,
        density: Density::Compact,
        accent: Accent::Teal,
        contrast: Contrast::High,
        colours: Colours::Safe,
    };
    assert_eq!(
        a.encode(ThemeMode::Dark),
        "theme dark\nscale 125\ndensity compact\naccent teal\ncontrast high\ncolours safe\n"
    );
}

#[test]
fn what_is_not_understood_leaves_that_setting_alone_and_spoils_no_other() {
    let text =
        "theme purple\nscale 5000\ndensity compact\naccent \ncontrast high extra\nwhatever 3\n\
                colours safe\nscale abc\nscale 150\nscale 99\n";
    let (theme, a) = Appearance::decode(text);
    assert_eq!(theme, ThemeMode::System, "a theme that does not exist");
    assert_eq!(
        a.scale, 150,
        "a size that is not one of the sizes, and not a number, are passed over"
    );
    assert_eq!(a.density, Density::Compact);
    assert_eq!(a.accent, Accent::Indigo, "a line with no value");
    assert_eq!(
        a.contrast,
        Contrast::Standard,
        "a line with more than a value"
    );
    assert_eq!(a.colours, Colours::Safe);
    assert!(is_default(Appearance::decode("").1));
    assert!(is_default(Appearance::decode("\0\u{1}garbage\n\n   \n").1));
}

#[test]
fn a_byte_order_mark_at_the_start_of_the_file_does_not_hide_the_first_line() {
    let (theme, a) = Appearance::decode("\u{feff}theme dark\nscale 125\n");
    assert_eq!((theme, a.scale), (ThemeMode::Dark, 125));
}

#[test]
fn settings_from_the_command_line_are_put_on_what_was_kept_and_leave_the_rest() {
    let kept = Appearance {
        accent: Accent::Teal,
        density: Density::Compact,
        ..Appearance::default()
    };
    let with = kept.with_spec("scale=150,contrast=high");
    assert_eq!((with.scale, with.contrast), (150, Contrast::High));
    assert_eq!(
        (with.accent, with.density),
        (Accent::Teal, Density::Compact)
    );
    assert_eq!(kept.with_spec(""), kept);
    assert_eq!(kept.with_spec("nonsense,scale=7"), kept);
}

#[test]
fn the_last_line_for_a_setting_counts() {
    assert_eq!(Appearance::decode("scale 90\nscale 175\n").1.scale, 175);
}

#[test]
fn the_command_line_names_settings_with_equals_signs_and_commas() {
    let a =
        Appearance::from_spec("scale=150,contrast=high,accent=green,colours=safe,density=compact");
    assert_eq!(
        a,
        Appearance {
            scale: 150,
            density: Density::Compact,
            accent: Accent::Green,
            contrast: Contrast::High,
            colours: Colours::Safe,
        }
    );
    assert!(is_default(Appearance::from_spec("")));
    assert!(is_default(Appearance::from_spec("scale=7,nonsense")));
}

#[test]
fn the_size_moves_one_step_at_a_time_and_stops_at_the_ends() {
    let mut a = Appearance::default();
    let mut seen = vec![a.scale];
    for _ in 0..10 {
        a = a.stepped(1);
        seen.push(a.scale);
    }
    assert_eq!(seen[..8], [100, 110, 125, 150, 175, 200, 200, 200]);
    assert_eq!(
        a.stepped(0).scale,
        100,
        "Cmd 0 goes back to the normal size"
    );
    for _ in 0..10 {
        a = a.stepped(-1);
    }
    assert_eq!(a.scale, 80, "and the smallest is the end at the other side");
    assert_eq!(a.stepped(-1).scale, 80);
    // the other settings are not touched
    let b = Appearance {
        accent: Accent::Blue,
        ..Appearance::default()
    };
    assert_eq!(b.stepped(1).accent, Accent::Blue);
    // a size that is not in the list (it cannot be read from a file, but could be set in code) steps from the normal one
    let odd = Appearance {
        scale: 133,
        ..Appearance::default()
    };
    assert_eq!(odd.stepped(1).scale, 110);
    assert_eq!(
        Appearance {
            scale: 125,
            ..Appearance::default()
        }
        .zoom(),
        1.25
    );
}

#[test]
fn every_choice_has_its_own_word_and_its_own_name() {
    fn distinct<T: Copy>(
        all: &[T],
        key: impl Fn(T) -> &'static str,
        label: impl Fn(T) -> &'static str,
    ) {
        let keys: std::collections::HashSet<_> = all.iter().map(|&c| key(c)).collect();
        let labels: std::collections::HashSet<_> = all.iter().map(|&c| label(c)).collect();
        assert_eq!((keys.len(), labels.len()), (all.len(), all.len()));
        assert!(all.iter().all(|&c| !key(c).contains(char::is_whitespace)));
    }
    distinct(Density::ALL, Density::key, Density::label);
    distinct(Accent::ALL, Accent::key, Accent::label);
    distinct(Contrast::ALL, Contrast::key, Contrast::label);
    distinct(Colours::ALL, Colours::key, Colours::label);
    distinct(&ThemeMode::ALL, ThemeMode::key, ThemeMode::label);
}

// ---- spacing ---------------------------------------------------------------------------------------------------------------------

#[test]
fn compact_is_tighter_but_no_target_is_smaller_than_the_least_the_guidelines_allow() {
    let (roomy, tight) = (Density::Comfortable.spacing(), Density::Compact.spacing());
    assert!(tight.item.y < roomy.item.y && tight.item.x < roomy.item.x);
    assert!(tight.button_padding.y < roomy.button_padding.y);
    assert!(tight.window_margin < roomy.window_margin);
    for spacing in [roomy, tight] {
        assert!(
            spacing.interact_height >= 24.0,
            "{}",
            spacing.interact_height
        );
    }
}

// ---- contrast ---------------------------------------------------------------------------------------------------------------------

#[test]
fn contrast_is_measured_as_the_guidelines_do() {
    let (black, white) = (Color32::BLACK, Color32::WHITE);
    assert!((contrast(black, white) - 21.0).abs() < 0.01);
    assert!(
        (contrast(white, black) - 21.0).abs() < 0.01,
        "either way round"
    );
    assert!((contrast(white, white) - 1.0).abs() < 1e-6);
    // the grey that is the lightest to pass level AA on white (4.5) is #767676, a known value of the guidelines' own examples
    let grey = |v| Color32::from_gray(v);
    assert!((contrast(grey(0x76), white) - 4.54).abs() < 0.02);
    assert!(contrast(grey(0x77), white) < 4.5);
}

/// Every colour that text or a line of the chrome is drawn on: the page, the cards, a box, a button under the pointer, and (last) a
/// chosen button's tint of the accent. The first four are the grounds of a line.
fn grounds(l: &Look) -> [(&'static str, Color32); 6] {
    [
        ("page", l.bg),
        ("surface", l.surface),
        ("alt surface", l.surface_alt),
        ("well", l.well),
        ("hover", l.hover),
        ("chosen", l.accent_soft),
    ]
}

/// `colour` as it looks over `behind` (a colour with some see-through in it, such as the accent shown through at a third).
fn flat(colour: Color32, behind: Color32) -> Color32 {
    let keep = 255 - u16::from(colour.a());
    let over = |c: u8, b: u8| (u16::from(c) + u16::from(b) * keep / 255).min(255) as u8;
    Color32::from_rgb(
        over(colour.r(), behind.r()),
        over(colour.g(), behind.g()),
        over(colour.b(), behind.b()),
    )
}

fn lum_pairs(l: &Look) -> [(&'static str, Color32, Color32); 4] {
    [
        ("text", l.text, l.surface),
        ("dim", l.dim, l.surface),
        ("accent", l.accent, l.surface),
        ("on_accent", l.on_accent, l.accent),
    ]
}

#[test]
fn every_accent_reads_on_the_page_and_carries_its_own_text_in_both_themes() {
    for &accent in Accent::ALL {
        for dark in [false, true] {
            let look = Appearance {
                accent,
                ..Appearance::default()
            }
            .look(dark);
            for (what, fg, bg) in lum_pairs(&look) {
                assert!(
                    contrast(fg, bg) >= READABLE,
                    "{what} of {accent:?}, dark={dark}: {:.2}",
                    contrast(fg, bg)
                );
            }
            // the tint that a chosen button or row is filled with keeps the ordinary text readable
            assert!(
                contrast(look.text, look.accent_soft) >= 7.0,
                "{accent:?}, dark={dark}"
            );
            assert_ne!(
                look.accent, look.accent2,
                "{accent:?}: the two colours are told apart"
            );
            assert!(!look.safe_colours);
        }
    }
}

#[test]
fn the_text_of_a_chosen_button_reads_on_its_ground_for_every_other_accent_and_for_high_contrast() {
    for &accent in Accent::ALL {
        for dark in [false, true] {
            for (mode, least) in [(Contrast::Standard, READABLE), (Contrast::High, 7.0)] {
                let look = Appearance {
                    accent,
                    contrast: mode,
                    ..Appearance::default()
                }
                .look(dark);
                if accent == Accent::Indigo && mode == Contrast::Standard {
                    // the designed look is left as it was: the accent shown through at a third, which is a paler ground than it
                    // should be (the text of a chosen button reads at about 3 to 1 on a dark page)
                    assert_eq!(look, design::look(dark));
                    continue;
                }
                // what the window draws a chosen button's text in is the colour of the selection's stroke
                let ground = flat(look.selection_fill, look.surface);
                assert!(
                    contrast(look.accent, ground) >= least,
                    "{accent:?} {mode:?} dark={dark}: {:.2}",
                    contrast(look.accent, ground)
                );
                assert_eq!(look.selection_fill.a(), 255, "an opaque ground");
            }
        }
    }
}

#[test]
fn the_accents_differ_from_each_other_and_only_the_indigo_is_the_designed_one() {
    for dark in [false, true] {
        let swatches: Vec<Color32> = Accent::ALL.iter().map(|a| a.swatch(dark)).collect();
        let unique: std::collections::HashSet<_> = swatches.iter().collect();
        assert_eq!(unique.len(), Accent::ALL.len(), "dark={dark}");
        assert_eq!(Accent::Indigo.swatch(dark), design::look(dark).accent);
    }
    // a swatch is what the chrome really uses
    let teal = Appearance {
        accent: Accent::Teal,
        ..Appearance::default()
    };
    assert_eq!(teal.look(false).accent, Accent::Teal.swatch(false));
}

#[test]
fn high_contrast_puts_text_lines_and_colours_far_enough_from_what_is_behind_them() {
    for &accent in Accent::ALL {
        for dark in [false, true] {
            let look = Appearance {
                accent,
                contrast: Contrast::High,
                ..Appearance::default()
            }
            .look(dark);
            for (what, colour) in [
                ("text", look.text),
                ("dim", look.dim),
                ("accent", look.accent),
                ("accent2", look.accent2),
                ("good", look.good),
                ("warn", look.warn),
                ("bad", look.bad),
            ] {
                for (page, bg) in grounds(&look) {
                    assert!(
                        contrast(colour, bg) >= 7.0,
                        "{what} on the {page}, {accent:?}, dark={dark}: {:.2}",
                        contrast(colour, bg)
                    );
                }
            }
            for (page, bg) in &grounds(&look)[..4] {
                assert!(
                    contrast(look.border, *bg) >= 3.0,
                    "the edge of a box on the {page}, {accent:?}, dark={dark}: {:.2}",
                    contrast(look.border, *bg)
                );
            }
            assert!(
                contrast(look.on_accent, look.accent) >= 7.0,
                "{accent:?}, dark={dark}"
            );
            assert!(look.outline > 1.0);
        }
    }
}

#[test]
fn the_standard_look_is_what_high_contrast_improves_on() {
    // so that the setting is not a change of nothing: the designed colours fall short of what the high-contrast look reaches
    for dark in [false, true] {
        let (standard, high) = (
            Appearance::default().look(dark),
            Appearance {
                contrast: Contrast::High,
                ..Appearance::default()
            }
            .look(dark),
        );
        assert!(
            contrast(standard.border, standard.surface) < 3.0,
            "dark={dark}"
        );
        assert!(contrast(standard.border, standard.bg) < 3.0, "dark={dark}");
        assert_ne!(standard, high);
        // never less than the designed colours had, and more where they fell short
        for (what, new, old) in [
            ("text", high.text, standard.text),
            ("dim", high.dim, standard.dim),
            ("accent", high.accent, standard.accent),
            ("border", high.border, standard.border),
        ] {
            assert!(
                contrast(new, high.surface) >= contrast(old, standard.surface),
                "{what}, dark={dark}"
            );
        }
        assert!(
            contrast(high.border, high.surface) > 2.0 * contrast(standard.border, standard.surface)
        );
    }
    // the dimmed text of the light theme is one that falls short of 7 until it is raised
    {
        let (standard, high) = (
            design::look(false),
            Appearance {
                contrast: Contrast::High,
                ..Appearance::default()
            }
            .look(false),
        );
        assert!(contrast(standard.dim, standard.surface) < 7.0);
        assert!(contrast(high.dim, high.surface) >= 7.0);
    }
}

#[test]
fn high_contrast_reaches_the_editor_too() {
    for dark in [false, true] {
        let plain = crate::ui::theme::palette(dark);
        let high = Appearance {
            contrast: Contrast::High,
            ..Appearance::default()
        }
        .palette(plain);
        assert!(contrast(high.dim, high.text_bg) >= 7.0, "dark={dark}");
        assert!(contrast(high.line, high.text_bg) >= 3.0, "dark={dark}");
        assert!(contrast(high.fg, high.text_bg) >= 7.0, "dark={dark}");
        assert_eq!(high.text_bg, plain.text_bg, "the page itself is as it was");
    }
}

#[test]
fn letters_are_always_in_a_colour_that_reads_on_their_cell() {
    for dark in [false, true] {
        for contrast_mode in [Contrast::Standard, Contrast::High] {
            for colours in [Colours::Standard, Colours::Safe] {
                let look = Appearance {
                    contrast: contrast_mode,
                    colours,
                    ..Appearance::default()
                }
                .look(dark);
                let fills = b"ACGTU".iter().filter_map(|&b| base_fill(b, &look)).chain(
                    b"ACDEFGHIKLMNPQRSTVWY"
                        .iter()
                        .filter_map(|&b| residue_fill(b, &look)),
                );
                for fill in fills {
                    let ink = design::ink_on(&look, fill);
                    assert!(
                        contrast(ink, fill) >= READABLE,
                        "{fill:?}, dark={dark}, {colours:?}: {:.2}",
                        contrast(ink, fill)
                    );
                }
            }
        }
    }
}

#[test]
fn letters_keep_the_text_colour_of_the_theme_wherever_that_reads() {
    let light = design::look(false);
    assert_eq!(design::ink_on(&light, Color32::WHITE), light.text);
    assert_eq!(
        design::ink_on(&light, Color32::from_rgb(0x9a, 0xd9, 0x9a)),
        light.text
    );
    // and give way to black or white where it does not
    assert_eq!(
        design::ink_on(&light, Color32::from_rgb(0x10, 0x10, 0x30)),
        Color32::WHITE
    );
    assert_eq!(
        design::ink_on(&design::look(true), Color32::from_rgb(0xf0, 0xf0, 0x80)),
        Color32::BLACK
    );
}

// ---- colour-blind-safe colours -----------------------------------------------------------------------------------------------------

/// How the colours look to people with each kind of colour blindness, and how far apart they are then: written out here with its own
/// matrices (Machado, Oliveira and Fernandes 2009, full strength, applied to linear sRGB) and CIE Lab (D65), so that it is not the code
/// that is under test that does the measuring. The numbers pinned below were produced by a separate script, in Python.
mod seen {
    use super::Color32;

    type Matrix = [[f64; 3]; 3];
    const NORMAL: Matrix = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    const PROTAN: Matrix = [
        [0.152286, 1.052583, -0.204868],
        [0.114503, 0.786281, 0.099216],
        [-0.003882, -0.048116, 1.051998],
    ];
    const DEUTAN: Matrix = [
        [0.367322, 0.860646, -0.227968],
        [0.280085, 0.672501, 0.047413],
        [-0.011820, 0.042940, 0.968881],
    ];
    const TRITAN: Matrix = [
        [1.255528, -0.076749, -0.178779],
        [-0.078411, 0.930809, 0.147602],
        [0.004733, 0.691367, 0.303900],
    ];
    pub const KINDS: [(&str, Matrix); 4] = [
        ("normal", NORMAL),
        ("protan", PROTAN),
        ("deutan", DEUTAN),
        ("tritan", TRITAN),
    ];
    const XYZ: Matrix = [
        [0.4124564, 0.3575761, 0.1804375],
        [0.2126729, 0.7151522, 0.0721750],
        [0.0193339, 0.1191920, 0.9503041],
    ];

    fn linear(v: u8) -> f64 {
        let v = f64::from(v) / 255.0;
        if v <= 0.04045 {
            v / 12.92
        } else {
            ((v + 0.055) / 1.055).powf(2.4)
        }
    }

    fn times(m: &Matrix, v: [f64; 3]) -> [f64; 3] {
        m.map(|row| row[0] * v[0] + row[1] * v[1] + row[2] * v[2])
    }

    pub fn lab(c: Color32, kind: &Matrix) -> [f64; 3] {
        let seen =
            times(kind, [linear(c.r()), linear(c.g()), linear(c.b())]).map(|v| v.clamp(0.0, 1.0));
        let white = times(&XYZ, [1.0, 1.0, 1.0]);
        let xyz = times(&XYZ, seen);
        let f = |t: f64| {
            if t > (6.0f64 / 29.0).powi(3) {
                t.cbrt()
            } else {
                t / (3.0 * (6.0f64 / 29.0).powi(2)) + 4.0 / 29.0
            }
        };
        let (x, y, z) = (
            f(xyz[0] / white[0]),
            f(xyz[1] / white[1]),
            f(xyz[2] / white[2]),
        );
        [116.0 * y - 16.0, 500.0 * (x - y), 200.0 * (y - z)]
    }

    /// The smallest distance between two of `colours`, for each kind of vision (in the order of `KINDS`).
    pub fn worst(colours: &[Color32]) -> [f64; 4] {
        KINDS.map(|(_, kind)| {
            let labs: Vec<_> = colours.iter().map(|&c| lab(c, &kind)).collect();
            let mut least = f64::INFINITY;
            for (i, a) in labs.iter().enumerate() {
                for b in &labs[i + 1..] {
                    least = least.min(
                        ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2))
                            .sqrt(),
                    );
                }
            }
            least
        })
    }
}

fn look_with(colours: Colours, dark: bool) -> Look {
    Appearance {
        colours,
        ..Appearance::default()
    }
    .look(dark)
}

/// The four kinds of colours that the viewers draw, for a theme and a choice of colours.
fn sets(colours: Colours, dark: bool) -> [(&'static str, Vec<Color32>); 3] {
    let look = look_with(colours, dark);
    // as the cells look on the page: the standard colours of the dark residues are see-through (a share of the colour, the rest the page)
    let on_page = |c: Color32| flat(c, look.surface);
    let bases = b"ACGT"
        .iter()
        .filter_map(|&b| base_fill(b, &look).map(on_page))
        .collect();
    let groups = b"LKDSCGPY"
        .iter()
        .filter_map(|&b| residue_fill(b, &look).map(on_page))
        .collect();
    let trace = b"ACGTN".iter().map(|&b| trace_color(b, &look)).collect();
    [("bases", bases), ("residues", groups), ("trace", trace)]
}

fn smallest(worst: [f64; 4]) -> f64 {
    worst.into_iter().fold(f64::INFINITY, f64::min)
}

#[test]
fn the_safe_colours_stay_apart_for_every_kind_of_colour_blindness_and_the_standard_ones_do_not() {
    // what was measured separately (see `seen`): the smallest distance for each of normal, protan, deutan and tritan vision
    let pinned: [(&str, bool, [f64; 4], [f64; 4]); 6] = [
        (
            "bases",
            false,
            [36.2, 8.9, 8.7, 11.6],
            [42.7, 32.4, 23.9, 14.2],
        ),
        (
            "bases",
            true,
            [42.3, 16.2, 13.9, 14.5],
            [61.8, 35.1, 30.6, 15.3],
        ),
        (
            "residues",
            false,
            [13.9, 8.0, 5.0, 4.2],
            [22.9, 19.3, 19.3, 19.3],
        ),
        (
            "residues",
            true,
            [16.6, 8.7, 4.6, 5.8],
            [16.3, 13.4, 13.4, 13.4],
        ),
        (
            "trace",
            false,
            [46.6, 4.3, 14.6, 30.9],
            [52.8, 22.7, 17.9, 17.0],
        ),
        (
            "trace",
            true,
            [46.6, 4.3, 14.6, 30.9],
            [41.6, 20.7, 17.9, 16.1],
        ),
    ];
    for (name, dark, standard, safe) in pinned {
        let find = |colours| {
            sets(colours, dark)
                .into_iter()
                .find(|(n, _)| *n == name)
                .map(|(_, c)| seen::worst(&c))
                .expect("a set")
        };
        let (std_now, safe_now) = (find(Colours::Standard), find(Colours::Safe));
        for (k, (kind, _)) in seen::KINDS.iter().enumerate() {
            assert!(
                (std_now[k] - standard[k]).abs() < 0.15,
                "{name} standard, {kind}, dark={dark}: {:.2} against {}",
                std_now[k],
                standard[k]
            );
            assert!(
                (safe_now[k] - safe[k]).abs() < 0.15,
                "{name} safe, {kind}, dark={dark}: {:.2} against {}",
                safe_now[k],
                safe[k]
            );
        }
        // the point of it: the worst case is much better, and never as low as the standard colours' worst
        assert!(
            smallest(safe_now) > smallest(std_now),
            "{name}, dark={dark}"
        );
        assert!(
            smallest(safe_now) >= 13.0,
            "{name}, dark={dark}: {:.1}",
            smallest(safe_now)
        );
    }
}

#[test]
fn the_safe_colours_are_a_choice_that_every_viewer_can_see_and_the_gap_still_has_none() {
    let safe = look_with(Colours::Safe, false);
    let plain = look_with(Colours::Standard, false);
    for dark in [false, true] {
        let (safe, plain) = (
            look_with(Colours::Safe, dark),
            look_with(Colours::Standard, dark),
        );
        for &b in b"ACGTUacgtu" {
            assert!(base_fill(b, &safe).is_some(), "{}", char::from(b));
            assert_ne!(
                base_fill(b, &safe),
                base_fill(b, &plain),
                "{}",
                char::from(b)
            );
        }
        for &aa in b"ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy" {
            assert!(residue_fill(aa, &safe).is_some(), "{}", char::from(aa));
        }
        for &none in b"-.*Xx?" {
            assert_eq!(residue_fill(none, &safe), None, "{}", char::from(none));
        }
        for &none in b"-.*Xx?" {
            assert_eq!(base_fill(none, &safe), None, "{}", char::from(none));
        }
        assert_eq!(
            base_fill(b'a', &safe),
            base_fill(b'A', &safe),
            "lower case reads like upper case"
        );
        assert_eq!(residue_fill(b'k', &safe), residue_fill(b'K', &safe));
        assert_eq!(
            base_fill(b'U', &safe),
            base_fill(b'T', &safe),
            "RNA reads like DNA"
        );
        assert_eq!(trace_color(b'a', &safe), trace_color(b'A', &safe));
    }
    // the four trace letters differ, and G stays dark or light by the theme as in the original
    let letters: std::collections::HashSet<_> =
        b"ACGT".iter().map(|&b| trace_color(b, &safe)).collect();
    assert_eq!(letters.len(), 4);
    assert_ne!(
        trace_color(b'G', &look_with(Colours::Safe, false)),
        trace_color(b'G', &look_with(Colours::Safe, true))
    );
    assert_eq!(trace_color(b'N', &safe), trace_color(b'-', &safe));
    assert_ne!(trace_color(b'A', &safe), trace_color(b'A', &plain));
}

#[test]
fn the_cross_check_is_right_about_what_it_measures() {
    // identical colours are no distance apart, black and white are the whole of the lightness scale, and a colour
    // blind person sees no difference between two colours on a line that the matrices keep (a grey has none to lose)
    assert!(seen::worst(&[Color32::RED, Color32::RED])[0] < 1e-9);
    let (black, white) = (
        seen::lab(Color32::BLACK, &seen::KINDS[0].1),
        seen::lab(Color32::WHITE, &seen::KINDS[0].1),
    );
    assert!(black[0].abs() < 1e-6 && (white[0] - 100.0).abs() < 1e-3);
    for (_, kind) in &seen::KINDS {
        let grey = seen::lab(Color32::from_gray(128), kind);
        assert!(grey[1].abs() < 0.5 && grey[2].abs() < 0.5, "{grey:?}");
    }
    // red and green are far apart for normal vision and close for a person who cannot tell them
    let both = [
        Color32::from_rgb(200, 40, 40),
        Color32::from_rgb(40, 150, 40),
    ];
    let worst = seen::worst(&both);
    assert!(worst[0] > 60.0 && worst[2] < 0.5 * worst[0], "{worst:?}");
}

#[test]
fn the_contig_view_and_the_editor_palette_take_the_safe_colours_too() {
    use crate::ui::windows::contig_view::base_color;
    for dark in [false, true] {
        let plain = Appearance::default().palette(crate::ui::theme::palette(dark));
        let safe = Appearance {
            colours: Colours::Safe,
            ..Appearance::default()
        }
        .palette(crate::ui::theme::palette(dark));
        assert!(!plain.safe_colours && safe.safe_colours);
        for &b in b"ACGT" {
            assert_ne!(
                base_color(b, &plain),
                base_color(b, &safe),
                "{}",
                char::from(b)
            );
            assert_eq!(base_color(b, &safe), safe::trace(b, dark));
        }
        // high contrast alone does not change them
        let high = Appearance {
            contrast: Contrast::High,
            ..Appearance::default()
        }
        .palette(crate::ui::theme::palette(dark));
        assert!(!high.safe_colours);
        assert_eq!(base_color(b'A', &high), base_color(b'A', &plain));
    }
}
