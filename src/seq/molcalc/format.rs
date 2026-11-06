//! Numbers in and out of the calculator the way the original's Tcl reads and writes them: a double is written with the
//! shortest digits that read back the same (`tcl_double`), `format %.3g` is C's (`g3`), and `round()` goes half away from zero.

/// A number as the fields accept it (`string is double`): digits with a sign, a point and an exponent, nothing else.
pub fn number(text: &str) -> Option<f64> {
    let t = text.trim();
    // Rust also reads "inf", "nan" and "infinity": a field never holds those
    if t.is_empty()
        || !t
            .bytes()
            .all(|b| b.is_ascii_digit() || matches!(b, b'+' | b'-' | b'.' | b'e' | b'E'))
    {
        return None;
    }
    t.parse::<f64>().ok().filter(|v| v.is_finite())
}

/// `round(100 * x) / 100.0`: Tcl's `round` goes half away from zero, and a result of -0 is 0 there (it is an integer).
pub fn round2(x: f64) -> f64 {
    ((100.0 * x).round() + 0.0) / 100.0
}

/// `format %.2f [expr {round(100 * x) / 100.0}]`.
pub fn fixed2(x: f64) -> String {
    format!("{:.2}", round2(x))
}

/// A double as Tcl prints it: the shortest digits that read back the same, `.0` after a whole number, and an exponent
/// only below 1e-4 or from 1e17 up (`1e-5`, `1e+17`).
pub fn tcl_double(v: f64) -> String {
    if v.is_nan() {
        return "NaN".into();
    }
    if v.is_infinite() {
        return if v > 0.0 { "Inf" } else { "-Inf" }.into();
    }
    let sign = if v.is_sign_negative() { "-" } else { "" };
    if v == 0.0 {
        return format!("{sign}0.0");
    }
    // the shortest digits come from the scientific form: "1.2345e4" is the digits 12345 and the exponent 4
    let sci = format!("{:e}", v.abs());
    let (mantissa, exp) = sci.split_once('e').unwrap_or((&sci, "0"));
    let exp: i32 = exp.parse().unwrap_or(0);
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    let body = if !(-4..17).contains(&exp) {
        let rest = &digits[1..];
        let point = if rest.is_empty() { "" } else { "." };
        format!(
            "{}{point}{rest}e{}{}",
            &digits[..1],
            if exp < 0 { "-" } else { "+" },
            exp.abs()
        )
    } else if exp >= 0 {
        let whole = exp as usize + 1;
        if digits.len() <= whole {
            format!("{digits}{}.0", "0".repeat(whole - digits.len()))
        } else {
            format!("{}.{}", &digits[..whole], &digits[whole..])
        }
    } else {
        format!("0.{}{digits}", "0".repeat((-exp - 1) as usize))
    };
    format!("{sign}{body}")
}

/// `format %.3g` as C writes it: three significant digits, no trailing zeros, an exponent below 1e-4 or from 1e3 up
/// (`1e+03`, `2.5e-07`), and `inf` and `nan` in lower case (Tcl's `format` hands the number to C's `printf`).
pub fn g3(v: f64) -> String {
    const PRECISION: i32 = 3;
    if v.is_nan() {
        return "nan".into();
    }
    if v.is_infinite() {
        return if v > 0.0 { "inf" } else { "-inf" }.into();
    }
    if v == 0.0 {
        return if v.is_sign_negative() { "-0" } else { "0" }.into();
    }
    // the exponent is the one of the number rounded to three digits (999.5 is 1e+03)
    let sci = format!("{:.*e}", (PRECISION - 1) as usize, v);
    let (mantissa, exp) = sci.split_once('e').unwrap_or((&sci, "0"));
    let exp: i32 = exp.parse().unwrap_or(0);
    if !(-4..PRECISION).contains(&exp) {
        format!(
            "{}e{}{:02}",
            trim_zeros(mantissa),
            if exp < 0 { "-" } else { "+" },
            exp.abs()
        )
    } else {
        let fixed = format!("{:.*}", (PRECISION - 1 - exp) as usize, v);
        trim_zeros(&fixed).to_string()
    }
}

/// The digits after a decimal point lose their trailing zeros, and the point with them if none are left.
fn trim_zeros(s: &str) -> &str {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        s
    }
}
