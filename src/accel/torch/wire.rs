//! The messages between the program and a PyTorch worker (`worker.py` describes them): an operation name and named arrays of
//! bytes, 32-bit and 64-bit numbers, little-endian, over the worker's standard input and output.

use std::io::{Read, Write};

/// The most bytes one array may hold when read back: a worker's answer is a matrix of counts or a number, never more than this.
const MOST_BYTES: u64 = 1 << 32;

#[derive(Clone, Debug, PartialEq)]
pub enum Array {
    U8(Vec<u8>),
    U32(Vec<u32>),
    F32(Vec<f32>),
    F64(Vec<f64>),
}

impl Array {
    fn code(&self) -> u8 {
        match self {
            Array::U8(_) => 1,
            Array::U32(_) => 2,
            Array::F32(_) => 3,
            Array::F64(_) => 4,
        }
    }

    fn len(&self) -> usize {
        match self {
            Array::U8(v) => v.len(),
            Array::U32(v) => v.len(),
            Array::F32(v) => v.len(),
            Array::F64(v) => v.len(),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Message {
    pub op: String,
    fields: Vec<(String, Array)>,
}

impl Message {
    pub fn new(op: &str) -> Message {
        Message {
            op: op.to_string(),
            fields: Vec::new(),
        }
    }

    /// This message with the array `name` added.
    pub fn with(mut self, name: &str, array: Array) -> Message {
        debug_assert!(name.len() < 256 && self.op.len() < 256);
        self.fields.push((name.to_string(), array));
        self
    }

    /// A whole number as an array of one.
    pub fn with_number(self, name: &str, value: usize) -> Message {
        self.with(name, Array::U32(vec![value as u32]))
    }

    pub fn array(&self, name: &str) -> Option<&Array> {
        self.fields.iter().find(|(n, _)| n == name).map(|(_, a)| a)
    }

    pub fn f64s(&self, name: &str) -> Option<&[f64]> {
        match self.array(name)? {
            Array::F64(v) => Some(v),
            _ => None,
        }
    }

    pub fn f32s(&self, name: &str) -> Option<&[f32]> {
        match self.array(name)? {
            Array::F32(v) => Some(v),
            _ => None,
        }
    }

    pub fn u8s(&self, name: &str) -> Option<&[u8]> {
        match self.array(name)? {
            Array::U8(v) => Some(v),
            _ => None,
        }
    }

    pub fn u32s(&self, name: &str) -> Option<&[u32]> {
        match self.array(name)? {
            Array::U32(v) => Some(v),
            _ => None,
        }
    }

    /// An array of bytes as text.
    pub fn text(&self, name: &str) -> Option<String> {
        match self.array(name)? {
            Array::U8(v) => Some(String::from_utf8_lossy(v).into_owned()),
            _ => None,
        }
    }

    pub fn write_to(&self, out: &mut impl Write) -> std::io::Result<()> {
        let mut bytes = vec![self.op.len() as u8];
        bytes.extend(self.op.as_bytes());
        bytes.extend((self.fields.len() as u16).to_le_bytes());
        for (name, array) in &self.fields {
            bytes.push(name.len() as u8);
            bytes.extend(name.as_bytes());
            bytes.push(array.code());
            bytes.extend((array.len() as u64).to_le_bytes());
            match array {
                Array::U8(v) => bytes.extend(v),
                Array::U32(v) => v.iter().for_each(|x| bytes.extend(x.to_le_bytes())),
                Array::F32(v) => v.iter().for_each(|x| bytes.extend(x.to_le_bytes())),
                Array::F64(v) => v.iter().for_each(|x| bytes.extend(x.to_le_bytes())),
            }
        }
        out.write_all(&bytes)?;
        out.flush()
    }

    pub fn read_from(input: &mut impl Read) -> std::io::Result<Message> {
        let op = read_name(input)?;
        let mut count = [0u8; 2];
        input.read_exact(&mut count)?;
        let mut message = Message::new(&op);
        for _ in 0..u16::from_le_bytes(count) {
            let name = read_name(input)?;
            let mut head = [0u8; 9];
            input.read_exact(&mut head)?;
            let n = u64::from_le_bytes(head[1..].try_into().expect("eight bytes"));
            let width = match head[0] {
                1 => 1,
                2 | 3 => 4,
                4 => 8,
                other => return Err(invalid(format!("an array of type {other}"))),
            };
            if n.saturating_mul(width) > MOST_BYTES {
                return Err(invalid(format!("an array of {n} numbers")));
            }
            let mut raw = vec![0u8; n as usize * width as usize];
            input.read_exact(&mut raw)?;
            let array = match head[0] {
                1 => Array::U8(raw),
                2 => Array::U32(numbers(&raw, u32::from_le_bytes)),
                3 => Array::F32(numbers(&raw, f32::from_le_bytes)),
                _ => Array::F64(numbers(&raw, f64::from_le_bytes)),
            };
            message = message.with(&name, array);
        }
        Ok(message)
    }
}

fn read_name(input: &mut impl Read) -> std::io::Result<String> {
    let mut length = [0u8; 1];
    input.read_exact(&mut length)?;
    let mut name = vec![0u8; length[0] as usize];
    input.read_exact(&mut name)?;
    String::from_utf8(name).map_err(|_| invalid("a name that is not text".into()))
}

fn numbers<T, const N: usize>(raw: &[u8], make: fn([u8; N]) -> T) -> Vec<T> {
    raw.as_chunks::<N>().0.iter().map(|c| make(*c)).collect()
}

fn invalid(what: String) -> std::io::Error {
    std::io::Error::new(
        std::io::ErrorKind::InvalidData,
        format!("The worker sent {what}."),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_message_reads_back_as_it_was_written() {
        let message = Message::new("counts")
            .with_number("first", 3)
            .with("weights", Array::U32(vec![0, 1, 1024, u32::MAX]))
            .with("f", Array::F32(vec![1.5, -0.0, f32::MIN_POSITIVE]))
            .with("d", Array::F64(vec![std::f64::consts::PI, -1e300]))
            .with("text", Array::U8(b"hello".to_vec()));
        let mut bytes = Vec::new();
        message.write_to(&mut bytes).expect("written");
        let back = Message::read_from(&mut bytes.as_slice()).expect("read");
        assert_eq!(back, message);
        assert_eq!(back.u32s("first"), Some([3].as_slice()));
        assert_eq!(back.text("text").as_deref(), Some("hello"));
        assert_eq!(back.f64s("missing"), None);
        // the wrong type for a name is not the right array
        assert_eq!(back.f64s("f"), None);
    }

    #[test]
    fn a_message_cut_short_or_with_nonsense_is_an_error_not_a_hang_or_a_panic() {
        let mut bytes = Vec::new();
        Message::new("x")
            .with("a", Array::F64(vec![1.0, 2.0]))
            .write_to(&mut bytes)
            .expect("written");
        for cut in 0..bytes.len() {
            assert!(Message::read_from(&mut &bytes[..cut]).is_err(), "{cut}");
        }
        // a type that does not exist, and an array that claims to be huge
        let mut odd = bytes.clone();
        odd[1 + 1 + 2 + 1 + 1] = 9;
        assert!(Message::read_from(&mut odd.as_slice()).is_err());
        let mut huge = bytes.clone();
        let at = 1 + 1 + 2 + 1 + 1 + 1;
        huge[at..at + 8].copy_from_slice(&u64::MAX.to_le_bytes());
        assert!(Message::read_from(&mut huge.as_slice()).is_err());
    }
}
