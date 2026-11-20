//! JSON, as far as the answers of GitHub need it: objects, arrays, strings (with every escape, `\u` pairs of surrogates too), numbers,
//! `true`, `false` and `null`. Nothing is taken that is nested more than `MOST_DEPTH` deep, or that has anything after the value.

const MOST_DEPTH: usize = 64;

#[derive(Clone, Debug, PartialEq)]
pub enum Json {
    Null,
    Bool(bool),
    Number(f64),
    String(String),
    Array(Vec<Json>),
    Object(Vec<(String, Json)>),
}

impl Json {
    pub fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Object(fields) => fields.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            Json::String(s) => Some(s),
            _ => None,
        }
    }

    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Json::Number(n) => Some(*n),
            _ => None,
        }
    }

    pub fn as_bool(&self) -> Option<bool> {
        match self {
            Json::Bool(b) => Some(*b),
            _ => None,
        }
    }

    pub fn as_array(&self) -> Option<&[Json]> {
        match self {
            Json::Array(items) => Some(items),
            _ => None,
        }
    }
}

struct Reader<'a> {
    text: &'a [u8],
    at: usize,
}

impl Reader<'_> {
    fn fail<T>(&self, what: &str) -> Result<T, String> {
        Err(format!("JSON: {what} at byte {}.", self.at))
    }

    fn skip_white(&mut self) {
        while self.at < self.text.len()
            && matches!(self.text[self.at], b' ' | b'\t' | b'\n' | b'\r')
        {
            self.at += 1;
        }
    }

    fn peek(&self) -> Option<u8> {
        self.text.get(self.at).copied()
    }

    fn word(&mut self, word: &str, value: Json) -> Result<Json, String> {
        if self.text[self.at..].starts_with(word.as_bytes()) {
            self.at += word.len();
            Ok(value)
        } else {
            self.fail("an unknown word")
        }
    }

    fn value(&mut self, depth: usize) -> Result<Json, String> {
        if depth > MOST_DEPTH {
            return self.fail("nested too deep");
        }
        self.skip_white();
        match self.peek() {
            None => self.fail("the text ends where a value should be"),
            Some(b'{') => self.object(depth),
            Some(b'[') => self.array(depth),
            Some(b'"') => Ok(Json::String(self.string()?)),
            Some(b't') => self.word("true", Json::Bool(true)),
            Some(b'f') => self.word("false", Json::Bool(false)),
            Some(b'n') => self.word("null", Json::Null),
            Some(b'-' | b'0'..=b'9') => self.number(),
            Some(_) => self.fail("a character that begins no value"),
        }
    }

    fn number(&mut self) -> Result<Json, String> {
        let start = self.at;
        while self.at < self.text.len()
            && matches!(
                self.text[self.at],
                b'-' | b'+' | b'.' | b'e' | b'E' | b'0'..=b'9'
            )
        {
            self.at += 1;
        }
        let slice = std::str::from_utf8(&self.text[start..self.at]).unwrap_or("");
        match slice.parse::<f64>() {
            Ok(n) if n.is_finite() => Ok(Json::Number(n)),
            _ => {
                self.at = start;
                self.fail("a number that is not one")
            }
        }
    }

    fn hex4(&mut self) -> Result<u32, String> {
        let digits = self
            .text
            .get(self.at..self.at + 4)
            .ok_or("JSON: a \\u escape is cut short.")?;
        let text = std::str::from_utf8(digits)
            .map_err(|_| "JSON: a \\u escape holds no hex digits.".to_string())?;
        let value = u32::from_str_radix(text, 16)
            .map_err(|_| "JSON: a \\u escape holds no hex digits.".to_string())?;
        self.at += 4;
        Ok(value)
    }

    fn string(&mut self) -> Result<String, String> {
        self.at += 1; // the opening quote
        let mut out: Vec<u8> = Vec::new();
        loop {
            let Some(c) = self.peek() else {
                return self.fail("a string that never ends");
            };
            self.at += 1;
            match c {
                b'"' => {
                    return String::from_utf8(out)
                        .map_err(|_| "JSON: a string is not UTF-8.".to_string())
                }
                b'\\' => {
                    let Some(e) = self.peek() else {
                        return self.fail("an escape that is cut short");
                    };
                    self.at += 1;
                    let mut push = |ch: char| out.extend_from_slice(ch.to_string().as_bytes());
                    match e {
                        b'"' => push('"'),
                        b'\\' => push('\\'),
                        b'/' => push('/'),
                        b'b' => push('\u{8}'),
                        b'f' => push('\u{c}'),
                        b'n' => push('\n'),
                        b'r' => push('\r'),
                        b't' => push('\t'),
                        b'u' => {
                            let high = self.hex4()?;
                            let code = if (0xD800..0xDC00).contains(&high) {
                                // the first half of a pair: the second half must follow
                                if self.text[self.at..].starts_with(b"\\u") {
                                    self.at += 2;
                                    let low = self.hex4()?;
                                    if !(0xDC00..0xE000).contains(&low) {
                                        return self
                                            .fail("half of a surrogate pair with no other half");
                                    }
                                    0x10000 + ((high - 0xD800) << 10) + (low - 0xDC00)
                                } else {
                                    return self
                                        .fail("half of a surrogate pair with no other half");
                                }
                            } else {
                                high
                            };
                            let Some(ch) = char::from_u32(code) else {
                                return self.fail("an escape that is no character");
                            };
                            out.extend_from_slice(ch.to_string().as_bytes());
                        }
                        _ => return self.fail("an escape that does not exist"),
                    }
                }
                c if c < 0x20 => return self.fail("a control character inside a string"),
                c => out.push(c),
            }
        }
    }

    fn array(&mut self, depth: usize) -> Result<Json, String> {
        self.at += 1;
        let mut items = Vec::new();
        self.skip_white();
        if self.peek() == Some(b']') {
            self.at += 1;
            return Ok(Json::Array(items));
        }
        loop {
            items.push(self.value(depth + 1)?);
            self.skip_white();
            match self.peek() {
                Some(b',') => self.at += 1,
                Some(b']') => {
                    self.at += 1;
                    return Ok(Json::Array(items));
                }
                _ => return self.fail("a comma or a closing bracket was expected"),
            }
        }
    }

    fn object(&mut self, depth: usize) -> Result<Json, String> {
        self.at += 1;
        let mut fields = Vec::new();
        self.skip_white();
        if self.peek() == Some(b'}') {
            self.at += 1;
            return Ok(Json::Object(fields));
        }
        loop {
            self.skip_white();
            if self.peek() != Some(b'"') {
                return self.fail("a name in quotes was expected");
            }
            let key = self.string()?;
            self.skip_white();
            if self.peek() != Some(b':') {
                return self.fail("a colon was expected");
            }
            self.at += 1;
            fields.push((key, self.value(depth + 1)?));
            self.skip_white();
            match self.peek() {
                Some(b',') => self.at += 1,
                Some(b'}') => {
                    self.at += 1;
                    return Ok(Json::Object(fields));
                }
                _ => return self.fail("a comma or a closing brace was expected"),
            }
        }
    }
}

/// The value of a JSON text; why not, when the text is not JSON or has anything after its value.
pub fn parse(text: &str) -> Result<Json, String> {
    let mut r = Reader {
        text: text.as_bytes(),
        at: 0,
    };
    let value = r.value(0)?;
    r.skip_white();
    if r.at != r.text.len() {
        return r.fail("something after the end of the value");
    }
    Ok(value)
}
