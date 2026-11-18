use super::*;
use std::collections::HashSet;

fn startup(view: &str) -> Startup {
    Startup {
        view: view.to_string(),
        ..Startup::default()
    }
}

fn few_frames(view: &str) -> ColorImage {
    draw(startup(view), 3, Duration::ZERO)
}

#[test]
fn the_program_is_drawn_without_a_window_at_the_size_of_the_window() {
    let image = few_frames("");
    assert_eq!(image.size, [2560, 1600]);
    let colors: HashSet<&Color32> = image.pixels.iter().collect();
    // menus, text, the map and the sequence: far more than a blank or half-drawn picture has
    assert!(colors.len() > 50, "only {} colours", colors.len());
}

#[test]
fn the_view_named_on_the_command_line_is_in_the_picture() {
    assert_ne!(few_frames("").pixels, few_frames("about").pixels);
}

#[test]
fn the_picture_is_saved_as_a_png_of_the_same_size() {
    let image = ColorImage::new([3, 2], vec![Color32::RED; 6]);
    let path = std::env::temp_dir().join(format!("helixloom-headless-{}.png", std::process::id()));
    let name = path.to_str().expect("a temporary path is text");
    save_png(name, &image).expect("the picture is saved");
    let file = std::io::BufReader::new(std::fs::File::open(&path).expect("the file is there"));
    let reader = png::Decoder::new(file).read_info().expect("it is a png");
    assert_eq!((reader.info().width, reader.info().height), (3, 2));
    std::fs::remove_file(&path).expect("the file is removed");
}

#[test]
fn the_window_and_the_fallback_never_both_take_the_picture() {
    assert!(claim_picture());
    assert!(!claim_picture());
}

fn zoomed(scale: u16) -> ColorImage {
    let start = Startup {
        appearance: Some(format!("scale={scale}")),
        ..Startup::default()
    };
    draw(start, 6, Duration::ZERO)
}

#[test]
fn a_zoomed_window_holds_fewer_points_and_the_picture_is_still_as_big_as_the_window() {
    let rect = frame_input_at(0, 2.0).screen_rect.expect("a screen");
    assert_eq!((rect.width(), rect.height()), (640.0, 400.0));
    let rect = frame_input_at(0, 0.8).screen_rect.expect("a screen");
    assert_eq!((rect.width(), rect.height()), (1600.0, 1000.0));
    let (large, small, normal) = (zoomed(200), zoomed(80), zoomed(100));
    for image in [&large, &small, &normal] {
        assert_eq!(image.size, [2560, 1600]);
        // the whole of it is drawn: a part that was left out would still be the black of an empty canvas
        let empty = image
            .pixels
            .iter()
            .filter(|c| **c == Color32::BLACK)
            .count();
        assert!(
            empty * 100 < image.pixels.len(),
            "{empty} of {} pixels are empty",
            image.pixels.len()
        );
    }
    assert_ne!(large.pixels, normal.pixels);
    assert_ne!(small.pixels, normal.pixels);
}
