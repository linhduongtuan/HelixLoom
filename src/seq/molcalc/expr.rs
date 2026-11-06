//! The Calculator of the original (Tools → Calculator, `calc` AppMain.tcl:34184): a line of arithmetic is made real by a
//! decimal point after every whole number, so that `1/3` is not 0 (`calc_update`, L34204, 0x18d6ad), and handed to Tcl's
//! `expr`. This does the same, with the part of `expr` that arithmetic needs: numbers (also `0x` hex and exponents),
//! `$pi` and `$e`, the operators of Tcl in its precedence (`**` is right to left and binds less tightly than a sign:
//! `-2**2` is 4), and the math functions. Whole numbers are 128-bit here; Tcl's go on to any size.
//!
//! What `expr` would run as a command or compare as text (`[…]`, quotes, `eq`, `ne`, `in`, `rand`) is not here: those
//! lines give the `------` that a line with a mistake in it gives.

use super::format::tcl_double;

/// What the result label says for a line that `expr` would refuse.
pub const MISTAKE: &str = "------";

/// The longest line that is worked out (in bytes) and the deepest nesting of brackets, calls, signs, powers and
/// conditions in it. A line beyond either is a mistake: nobody types one, and without the limits a long line of brackets
/// pasted into the window would use up the stack. (Not in the original, whose `expr` has limits of its own.)
pub const MAX_LINE: usize = 1000;
const MAX_DEPTH: usize = 100;

/// The result of a line of arithmetic, as the original's label shows it: `.0` is dropped from a whole real (`10.0` is
/// `10`), and a mistake gives `------`.
pub fn calculate(text: &str) -> String {
    if text.len() > MAX_LINE {
        return MISTAKE.into();
    }
    match evaluate(&make_real(text)) {
        Some(value) => {
            let shown = match value {
                Value::Int(i) => i.to_string(),
                Value::Real(r) => tcl_double(r),
            };
            shown
                .strip_suffix(".0")
                .map_or(shown.clone(), str::to_string)
        }
        None => MISTAKE.into(),
    }
}

/// `calc_update`: a point goes after every run of digits, unless the digits are part of a real number (after a point, before
/// one, in an exponent), a hex number, or next to a bit operator (`<<`, `>>`, `^`, `|`, `&`), where whole numbers are wanted.
pub fn make_real(text: &str) -> String {
    let chars: Vec<char> = text.chars().collect();
    let at =
        |i: isize| -> Option<char> { usize::try_from(i).ok().and_then(|i| chars.get(i)).copied() };
    // `string range`: the part from first to last, with both ends held to the string
    let range = |first: isize, last: isize| -> String {
        let (first, last) = (first.max(0), last.min(chars.len() as isize - 1));
        if first > last {
            return String::new();
        }
        chars[first as usize..=last as usize].iter().collect()
    };
    let upper = |s: String| s.to_uppercase();
    let mut points = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        if !chars[i].is_ascii_digit() {
            i += 1;
            continue;
        }
        let (a, mut b) = (i as isize, i as isize);
        while at(b + 1).is_some_and(|c| c.is_ascii_digit()) {
            b += 1;
        }
        i = b as usize + 1;
        let before = range(a - 2, a - 1);
        let after2 = range(b + 1, b + 2);
        let (before1, after1) = (range(a - 1, a - 1), range(b + 1, b + 1));
        let wanted = at(a - 1) != Some('.')
            && at(b + 1) != Some('.')
            && upper(before1.clone()) != "E"
            && !["E-", "E+", "0X"].contains(&upper(before.clone()).as_str())
            && upper(after1.clone()) != "X"
            && !["<<", ">>"].contains(&before.as_str())
            && !["^", "|", "&"].contains(&before1.as_str())
            && !["<<", ">>"].contains(&after2.as_str())
            && !["^", "|", "&"].contains(&after1.as_str());
        if wanted {
            points.push(b as usize);
        }
    }
    let mut out = chars;
    for (inserted, b) in points.into_iter().enumerate() {
        out.insert(inserted + b + 1, '.');
    }
    out.into_iter().collect()
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Value {
    Int(i128),
    Real(f64),
}

impl Value {
    fn real(self) -> f64 {
        match self {
            Value::Int(i) => i as f64,
            Value::Real(r) => r,
        }
    }

    fn truthy(self) -> bool {
        match self {
            Value::Int(i) => i != 0,
            Value::Real(r) => r != 0.0,
        }
    }
}

fn boolean(b: bool) -> Value {
    Value::Int(i128::from(b))
}

/// A real result: a number that is not a number (0/0, the square root of -1) is a mistake, infinity is not.
fn finite_or_inf(r: f64) -> Option<Value> {
    (!r.is_nan()).then_some(Value::Real(r))
}

#[derive(Clone, Debug, PartialEq)]
enum Token {
    Number(Value),
    Name(String),
    Variable(String),
    Op(&'static str),
}

const OPERATORS: [&str; 26] = [
    "**", "<<", ">>", "<=", ">=", "==", "!=", "&&", "||", "+", "-", "*", "/", "%", "(", ")", ",",
    "<", ">", "&", "|", "^", "~", "!", "?", ":",
];

/// The words of a line, or none if it has a character that `expr` would not take.
fn tokens(text: &str) -> Option<Vec<Token>> {
    let chars: Vec<char> = text.chars().collect();
    let mut out = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        if c.is_whitespace() {
            i += 1;
        } else if c.is_ascii_digit() || c == '.' {
            let (value, next) = number(&chars, i)?;
            out.push(Token::Number(value));
            i = next;
        } else if c.is_ascii_alphabetic() || c == '_' || c == '$' {
            let start = i + usize::from(c == '$');
            let mut end = start;
            while chars
                .get(end)
                .is_some_and(|c| c.is_ascii_alphanumeric() || *c == '_')
            {
                end += 1;
            }
            let word: String = chars[start..end].iter().collect();
            if word.is_empty() {
                return None;
            }
            out.push(if c == '$' {
                Token::Variable(word)
            } else {
                Token::Name(word)
            });
            i = end;
        } else {
            let rest: String = chars[i..chars.len().min(i + 2)].iter().collect();
            let op = *OPERATORS.iter().find(|op| rest.starts_with(**op))?;
            out.push(Token::Op(op));
            i += op.len();
        }
    }
    Some(out)
}

/// A number as Tcl reads one: `0x1F`, `12`, `1.`, `.5`, `1.5e-3`; and the place after it.
fn number(chars: &[char], start: usize) -> Option<(Value, usize)> {
    let digit = |i: usize| chars.get(i).is_some_and(char::is_ascii_digit);
    if chars[start] == '0' && matches!(chars.get(start + 1), Some('x' | 'X')) {
        let mut end = start + 2;
        while chars.get(end).is_some_and(char::is_ascii_hexdigit) {
            end += 1;
        }
        let text: String = chars[start + 2..end].iter().collect();
        return i128::from_str_radix(&text, 16)
            .ok()
            .map(|v| (Value::Int(v), end));
    }
    let mut end = start;
    while digit(end) {
        end += 1;
    }
    let mut real = false;
    if chars.get(end) == Some(&'.') {
        real = true;
        end += 1;
        while digit(end) {
            end += 1;
        }
    }
    if end == start || (real && end == start + 1) {
        return None;
    }
    // an exponent only counts with digits after it: `1e` is the number 1 and then a word
    if matches!(chars.get(end), Some('e' | 'E')) {
        let sign = usize::from(matches!(chars.get(end + 1), Some('+' | '-')));
        if digit(end + 1 + sign) {
            real = true;
            end += 1 + sign;
            while digit(end) {
                end += 1;
            }
        }
    }
    let text: String = chars[start..end].iter().collect();
    let value = if real {
        Value::Real(text.parse().ok()?)
    } else {
        Value::Int(text.parse().ok()?)
    };
    Some((value, end))
}

#[derive(Debug)]
enum Expr {
    Number(Value),
    Unary(&'static str, Box<Expr>),
    Binary(&'static str, Box<Expr>, Box<Expr>),
    Condition(Box<Expr>, Box<Expr>, Box<Expr>),
    Call(String, Vec<Expr>),
}

struct Parser {
    tokens: Vec<Token>,
    at: usize,
    /// How many nested parts (brackets, call arguments, signs, powers, conditions) are open.
    depth: usize,
}

/// The levels of binary operators from the loosest to the tightest (`**` and the unary signs have their own): `||`, `&&`,
/// `|`, `^`, `&`, equality, comparison, shifts, sums, products.
const LEVELS: [&[&str]; 9] = [
    &["||"],
    &["&&"],
    &["|"],
    &["^"],
    &["&"],
    &["==", "!="],
    &["<", ">", "<=", ">="],
    &["<<", ">>"],
    &["+", "-"],
];

impl Parser {
    fn peek(&self) -> Option<&Token> {
        self.tokens.get(self.at)
    }

    fn eat(&mut self, op: &str) -> bool {
        let found = matches!(self.peek(), Some(Token::Op(o)) if *o == op);
        if found {
            self.at += 1;
        }
        found
    }

    fn operator(&mut self, ops: &[&'static str]) -> Option<&'static str> {
        let Some(Token::Op(found)) = self.peek() else {
            return None;
        };
        let op = *ops.iter().find(|o| **o == *found)?;
        self.at += 1;
        Some(op)
    }

    /// Reads one more level of nesting with `parse`; too deep is a mistake.
    fn nested<T>(&mut self, parse: fn(&mut Self) -> Option<T>) -> Option<T> {
        if self.depth == MAX_DEPTH {
            return None;
        }
        self.depth += 1;
        let read = parse(self);
        self.depth -= 1;
        read
    }

    /// `x ? y : z`, the loosest, from right to left.
    fn condition(&mut self) -> Option<Expr> {
        let test = self.level(0)?;
        if !self.eat("?") {
            return Some(test);
        }
        let yes = self.nested(Self::condition)?;
        if !self.eat(":") {
            return None;
        }
        let no = self.nested(Self::condition)?;
        Some(Expr::Condition(Box::new(test), Box::new(yes), Box::new(no)))
    }

    /// One level of the binary operators, from left to right; the last one is the products.
    fn level(&mut self, n: usize) -> Option<Expr> {
        if n == LEVELS.len() {
            return self.product();
        }
        let mut left = self.level(n + 1)?;
        while let Some(op) = self.operator(LEVELS[n]) {
            left = Expr::Binary(op, Box::new(left), Box::new(self.level(n + 1)?));
        }
        Some(left)
    }

    fn product(&mut self) -> Option<Expr> {
        let mut left = self.power()?;
        while let Some(op) = self.operator(&["*", "/", "%"]) {
            left = Expr::Binary(op, Box::new(left), Box::new(self.power()?));
        }
        Some(left)
    }

    /// `**` goes from right to left, and a sign binds more tightly than it.
    fn power(&mut self) -> Option<Expr> {
        let base = self.unary()?;
        match self.operator(&["**"]) {
            Some(op) => {
                let exponent = self.nested(Self::power)?;
                Some(Expr::Binary(op, Box::new(base), Box::new(exponent)))
            }
            None => Some(base),
        }
    }

    fn unary(&mut self) -> Option<Expr> {
        match self.operator(&["-", "+", "~", "!"]) {
            Some(op) => Some(Expr::Unary(op, Box::new(self.nested(Self::unary)?))),
            None => self.primary(),
        }
    }

    fn primary(&mut self) -> Option<Expr> {
        let token = self.tokens.get(self.at)?.clone();
        self.at += 1;
        match token {
            Token::Number(v) => Some(Expr::Number(v)),
            Token::Variable(name) => match name.as_str() {
                "pi" => Some(Expr::Number(Value::Real(std::f64::consts::PI))),
                "e" => Some(Expr::Number(Value::Real(std::f64::consts::E))),
                _ => None,
            },
            Token::Op("(") => {
                let inner = self.nested(Self::condition)?;
                self.eat(")").then_some(inner)
            }
            Token::Name(name) => {
                if !self.eat("(") {
                    return None;
                }
                let mut args = Vec::new();
                if !self.eat(")") {
                    loop {
                        args.push(self.nested(Self::condition)?);
                        if self.eat(")") {
                            break;
                        }
                        if !self.eat(",") {
                            return None;
                        }
                    }
                }
                Some(Expr::Call(name, args))
            }
            Token::Op(_) => None,
        }
    }
}

/// What `expr` makes of a line; none for a mistake.
fn evaluate(text: &str) -> Option<Value> {
    let mut parser = Parser {
        tokens: tokens(text)?,
        at: 0,
        depth: 0,
    };
    let expr = parser.condition()?;
    (parser.at == parser.tokens.len()).then_some(())?;
    value(&expr)
}

fn value(expr: &Expr) -> Option<Value> {
    match expr {
        Expr::Number(v) => Some(*v),
        Expr::Unary(op, inner) => unary(op, value(inner)?),
        Expr::Binary("&&", a, b) => Some(boolean(value(a)?.truthy() && value(b)?.truthy())),
        Expr::Binary("||", a, b) => Some(boolean(value(a)?.truthy() || value(b)?.truthy())),
        Expr::Binary(op, a, b) => binary(op, value(a)?, value(b)?),
        Expr::Condition(test, yes, no) => {
            if value(test)?.truthy() {
                value(yes)
            } else {
                value(no)
            }
        }
        Expr::Call(name, args) => {
            let args: Vec<Value> = args.iter().map(value).collect::<Option<_>>()?;
            call(name, &args)
        }
    }
}

fn unary(op: &str, v: Value) -> Option<Value> {
    match (op, v) {
        ("-", Value::Int(i)) => i.checked_neg().map(Value::Int),
        ("-", Value::Real(r)) => Some(Value::Real(-r)),
        ("+", v) => Some(v),
        ("~", Value::Int(i)) => Some(Value::Int(!i)),
        ("!", v) => Some(boolean(!v.truthy())),
        _ => None,
    }
}

fn binary(op: &str, a: Value, b: Value) -> Option<Value> {
    // whole numbers stay whole under + - * / % ** and the bit operators; one real makes the result real
    if let (Value::Int(x), Value::Int(y)) = (a, b) {
        return whole(op, x, y);
    }
    let (x, y) = (a.real(), b.real());
    match op {
        "+" => finite_or_inf(x + y),
        "-" => finite_or_inf(x - y),
        "*" => finite_or_inf(x * y),
        "/" => finite_or_inf(x / y),
        // zero to a negative power is an error of the operator (`pow(0,-1)`, the function, is infinity)
        "**" if x == 0.0 && y < 0.0 => None,
        "**" => finite_or_inf(x.powf(y)),
        "<" => Some(boolean(x < y)),
        ">" => Some(boolean(x > y)),
        "<=" => Some(boolean(x <= y)),
        ">=" => Some(boolean(x >= y)),
        "==" => Some(boolean(x == y)),
        "!=" => Some(boolean(x != y)),
        // % << >> & ^ | are for whole numbers only
        _ => None,
    }
}

fn whole(op: &str, x: i128, y: i128) -> Option<Value> {
    let int = |v: Option<i128>| v.map(Value::Int);
    // a whole number that grows past 128 bits (Tcl's go on) goes on as a real number
    let grown = || binary(op, Value::Real(x as f64), Value::Real(y as f64));
    match op {
        "+" => int(x.checked_add(y)).or_else(grown),
        "-" => int(x.checked_sub(y)).or_else(grown),
        "*" => int(x.checked_mul(y)).or_else(grown),
        // Tcl rounds a quotient and a remainder down: -7/2 is -4, and -7 % 3 is 2
        "/" => int(floor_div(x, y)),
        "%" => int(floor_div(x, y).and_then(|q| x.checked_sub(y.checked_mul(q)?))),
        "**" if x == 0 && y < 0 => None,
        "**" => power(x, y).or_else(grown),
        // nothing shifted is still nothing
        "<<" if x == 0 && y >= 0 => Some(Value::Int(0)),
        "<<" if y < 0 => None,
        "<<" => int(u32::try_from(y)
            .ok()
            .filter(|n| *n < 127)
            .and_then(|n| x.checked_mul(1i128 << n)))
        .or_else(|| finite_or_inf(x as f64 * 2f64.powf(y as f64))),
        ">>" => u32::try_from(y).ok().map(|n| Value::Int(x >> n.min(127))),
        "&" => Some(Value::Int(x & y)),
        "^" => Some(Value::Int(x ^ y)),
        "|" => Some(Value::Int(x | y)),
        "<" => Some(boolean(x < y)),
        ">" => Some(boolean(x > y)),
        "<=" => Some(boolean(x <= y)),
        ">=" => Some(boolean(x >= y)),
        "==" => Some(boolean(x == y)),
        "!=" => Some(boolean(x != y)),
        _ => None,
    }
}

/// A quotient rounded down; none for a division by nothing.
fn floor_div(x: i128, y: i128) -> Option<i128> {
    let q = x.checked_div(y)?;
    Some(if x.checked_rem(y)? != 0 && ((x < 0) != (y < 0)) {
        q - 1
    } else {
        q
    })
}

/// A whole number to a whole power: a negative power is a fraction, so it is 0 (or ±1 for a base of ±1).
fn power(x: i128, y: i128) -> Option<Value> {
    if y >= 0 {
        return u32::try_from(y)
            .ok()
            .and_then(|n| x.checked_pow(n))
            .map(Value::Int);
    }
    match x {
        0 => None,
        1 => Some(Value::Int(1)),
        -1 => Some(Value::Int(if y % 2 == 0 { 1 } else { -1 })),
        _ => Some(Value::Int(0)),
    }
}

/// The math functions of `expr`: those that care whether a number is whole, and those that work on real numbers.
fn call(name: &str, args: &[Value]) -> Option<Value> {
    match name {
        "abs" | "int" | "entier" | "wide" | "round" | "isqrt" | "bool" | "max" | "min" => {
            by_kind(name, args)
        }
        _ => by_real(name, args),
    }
}

/// The functions that work on real numbers (a whole number is made real first).
fn by_real(name: &str, args: &[Value]) -> Option<Value> {
    let one = |f: fn(f64) -> f64| match args {
        [a] => finite_or_inf(f(a.real())),
        _ => None,
    };
    let two = |f: fn(f64, f64) -> f64| match args {
        [a, b] => finite_or_inf(f(a.real(), b.real())),
        _ => None,
    };
    let in_range = |lo: f64, hi: f64| matches!(args, [a] if (lo..=hi).contains(&a.real()));
    match name {
        "acos" => in_range(-1.0, 1.0).then(|| one(f64::acos)).flatten(),
        "asin" => in_range(-1.0, 1.0).then(|| one(f64::asin)).flatten(),
        "atan" => one(f64::atan),
        "atan2" => two(f64::atan2),
        "ceil" => one(f64::ceil),
        "cos" => one(f64::cos),
        "cosh" => one(f64::cosh),
        "double" => one(|x| x),
        "exp" => one(f64::exp),
        "floor" => one(f64::floor),
        "fmod" => match args {
            [_, b] if b.real() == 0.0 => None,
            _ => two(|x, y| x % y),
        },
        "hypot" => two(f64::hypot),
        "log" => in_range(0.0, f64::INFINITY).then(|| one(f64::ln)).flatten(),
        "log10" => in_range(0.0, f64::INFINITY)
            .then(|| one(f64::log10))
            .flatten(),
        "pow" => two(f64::powf),
        "sin" => one(f64::sin),
        "sinh" => one(f64::sinh),
        "sqrt" => in_range(0.0, f64::INFINITY)
            .then(|| one(f64::sqrt))
            .flatten(),
        "tan" => one(f64::tan),
        "tanh" => one(f64::tanh),
        _ => None,
    }
}

/// The functions that give a different answer for a whole number and for a real one.
fn by_kind(name: &str, args: &[Value]) -> Option<Value> {
    match name {
        "abs" => match args {
            [Value::Int(i)] => i.checked_abs().map(Value::Int),
            [Value::Real(r)] => Some(Value::Real(r.abs())),
            _ => None,
        },
        "int" | "entier" | "wide" => match args {
            [Value::Int(i)] => Some(Value::Int(*i)),
            [Value::Real(r)] if r.is_finite() && r.abs() < 1e38 => {
                Some(Value::Int(r.trunc() as i128))
            }
            _ => None,
        },
        "round" => match args {
            [Value::Int(i)] => Some(Value::Int(*i)),
            [Value::Real(r)] if r.is_finite() && r.abs() < 1e38 => {
                Some(Value::Int(r.round() as i128))
            }
            _ => None,
        },
        "isqrt" => match args {
            [Value::Int(i)] if *i >= 0 => Some(Value::Int((*i as f64).sqrt().floor() as i128)),
            [Value::Real(r)] if *r >= 0.0 && r.is_finite() => {
                Some(Value::Int(r.sqrt().floor() as i128))
            }
            _ => None,
        },
        "bool" => match args {
            [a] => Some(boolean(a.truthy())),
            _ => None,
        },
        "max" | "min" => {
            let (first, rest) = args.split_first()?;
            let better = |a: &Value, b: &Value| {
                if name == "max" {
                    b.real() > a.real()
                } else {
                    b.real() < a.real()
                }
            };
            Some(
                rest.iter()
                    .fold(*first, |best, v| if better(&best, v) { *v } else { best }),
            )
        }
        _ => None,
    }
}
