//! Sanger traces: the type that holds one (the four colour traces of a run, with the bases the machine called, where
//! each sits on the trace and how sure it was) and the reader of ABIF files (`.ab1`), the format of the Applied
//! Biosystems sequencers. SCF files (`scf.rs`) read into the same type, and `read` tells the two apart.
//!
//! The file is a header, then blocks of data, then a directory of 28-byte entries that say where each block is
//! (all numbers big-endian). Written from the published format description (Applied Biosystems, "ABIF File Format
//! Specification") and from the original's own reader, `read_abi_info` (AppMain.tcl:30572, 0x1684b9), which it
//! agrees with on real files from a 310, a 3100 and a 3730; checked against Biopython's reader too.

const MAGIC: &[u8; 4] = b"ABIF";
const HEADER: usize = 6;
const ENTRY: usize = 28;
/// Entries whose data is no longer than this sit in the entry itself, in place of the data's offset.
const INLINE: usize = 4;

/// The kinds of data an entry can hold (the element types of the specification) that this reader uses.
const BYTE: i16 = 1;
const CHAR: i16 = 2;
const WORD: i16 = 3;
const SHORT: i16 = 4;
const PSTRING: i16 = 18;
const CSTRING: i16 = 19;

/// Channels 9 to 12 hold the analysed signal, in the order of `FWO_`; 1 to 4 are the raw one.
const ANALYSED: [i32; 4] = [9, 10, 11, 12];
const RAW: [i32; 4] = [1, 2, 3, 4];

#[derive(Clone, Copy)]
struct Entry {
    name: [u8; 4],
    number: i32,
    kind: i16,
    size: i32,
    /// Where the data is, or (for data of up to four bytes) the data itself.
    offset: [u8; 4],
}

struct Abif<'a> {
    bytes: &'a [u8],
    entries: Vec<Entry>,
}

fn be_i16(b: &[u8]) -> i16 {
    i16::from_be_bytes([b[0], b[1]])
}

fn be_i32(b: &[u8]) -> i32 {
    i32::from_be_bytes([b[0], b[1], b[2], b[3]])
}

impl<'a> Abif<'a> {
    fn open(bytes: &'a [u8]) -> Result<Abif<'a>, String> {
        if bytes.len() < HEADER + ENTRY || &bytes[..4] != MAGIC {
            return Err(
                "This is not an ABIF trace file (.ab1): it does not start with \"ABIF\".".into(),
            );
        }
        // the header holds the entry that describes the directory
        let (count, at) = (be_i32(&bytes[HEADER + 12..]), be_i32(&bytes[HEADER + 20..]));
        let (count, at) = (usize::try_from(count), usize::try_from(at));
        let (Ok(count), Ok(at)) = (count, at) else {
            return Err("The trace file's directory is damaged.".into());
        };
        let end = count
            .checked_mul(ENTRY)
            .and_then(|n| n.checked_add(at))
            .filter(|&end| end <= bytes.len())
            .ok_or("The trace file is cut short: its directory runs past the end.")?;
        let entries = bytes[at..end]
            .as_chunks::<ENTRY>()
            .0
            .iter()
            .map(|e| Entry {
                name: [e[0], e[1], e[2], e[3]],
                number: be_i32(&e[4..]),
                kind: be_i16(&e[8..]),
                size: be_i32(&e[16..]),
                offset: [e[20], e[21], e[22], e[23]],
            })
            .collect();
        Ok(Abif { bytes, entries })
    }

    /// The data of entry `name` number `number`, copied.
    fn data(&self, name: &[u8; 4], number: i32) -> Option<Vec<u8>> {
        let e = self
            .entries
            .iter()
            .find(|e| &e.name == name && e.number == number)?;
        let size = usize::try_from(e.size).ok()?;
        if size <= INLINE {
            return Some(e.offset[..size].to_vec());
        }
        let at = usize::try_from(be_i32(&e.offset)).ok()?;
        self.bytes
            .get(at..at.checked_add(size)?)
            .map(<[u8]>::to_vec)
    }

    fn kind(&self, name: &[u8; 4], number: i32) -> Option<i16> {
        self.entries
            .iter()
            .find(|e| &e.name == name && e.number == number)
            .map(|e| e.kind)
    }

    fn shorts(&self, name: &[u8; 4], number: i32) -> Option<Vec<i16>> {
        let data = self.data(name, number)?;
        matches!(self.kind(name, number), Some(SHORT | WORD)).then(|| {
            let (pairs, _) = data.as_chunks::<2>();
            pairs.iter().map(|b| i16::from_be_bytes(*b)).collect()
        })
    }

    fn words(&self, name: &[u8; 4], number: i32) -> Option<Vec<u16>> {
        let data = self.data(name, number)?;
        matches!(self.kind(name, number), Some(SHORT | WORD)).then(|| {
            let (pairs, _) = data.as_chunks::<2>();
            pairs.iter().map(|b| u16::from_be_bytes(*b)).collect()
        })
    }

    fn bytes(&self, name: &[u8; 4], number: i32) -> Option<Vec<u8>> {
        matches!(self.kind(name, number), Some(BYTE | CHAR))
            .then(|| self.data(name, number))
            .flatten()
    }

    /// Text: a Pascal string (a length byte first), a C string (ends at a zero) or a plain array of characters.
    fn text(&self, name: &[u8; 4], number: i32) -> Option<String> {
        let data = self.data(name, number)?;
        let shown = match self.kind(name, number)? {
            PSTRING => {
                let n = usize::from(*data.first()?);
                data.get(1..=n.min(data.len().saturating_sub(1)))?.to_vec()
            }
            CSTRING => data.split(|&c| c == 0).next()?.to_vec(),
            CHAR => data.iter().copied().take_while(|&c| c != 0).collect(),
            _ => return None,
        };
        Some(String::from_utf8_lossy(&shown).trim().to_string())
    }
}

/// A moment a run began or ended, as the file gives it (`RUND` and `RUNT`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Stamp {
    pub year: u16,
    pub month: u8,
    pub day: u8,
    pub hour: u8,
    pub minute: u8,
    pub second: u8,
}

/// What the original shows for a run that gives no date or time: "01/01/00" and "00:00:00".
impl Default for Stamp {
    fn default() -> Self {
        Stamp {
            year: 2000,
            month: 1,
            day: 1,
            hour: 0,
            minute: 0,
            second: 0,
        }
    }
}

const MONTHS: [&str; 12] = [
    "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
];
const WEEKDAYS: [&str; 7] = ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"];

impl Stamp {
    /// The date as the original writes it: month, day and year, no padding ("2/19/2009").
    pub fn date_text(&self) -> String {
        format!("{}/{}/{}", self.month, self.day, self.year)
    }

    /// The time as the original writes it: the hour as it is, minutes and seconds with two digits ("1:19:30").
    pub fn time_text(&self) -> String {
        format!("{}:{:02}:{:02}", self.hour, self.minute, self.second)
    }

    /// The day of the week, 0 for Sunday (Sakamoto's method).
    fn weekday(&self) -> usize {
        const OFFSET: [usize; 12] = [0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4];
        let month = usize::from(self.month.clamp(1, 12));
        let year = usize::from(self.year) - usize::from(month < 3);
        (year + year / 4 - year / 100 + year / 400 + OFFSET[month - 1] + usize::from(self.day)) % 7
    }

    /// The moment in words, as the original's File Info shows it (its format is `%l:%M%p %a,%b %d, %Y`; a space
    /// is added after the comma): "1:19AM Thu, Feb 19, 2009".
    pub fn long_text(&self) -> String {
        let hour = match self.hour % 12 {
            0 => 12,
            h => h,
        };
        format!(
            "{hour}:{:02}{} {}, {} {:02}, {}",
            self.minute,
            if self.hour < 12 { "AM" } else { "PM" },
            WEEKDAYS[self.weekday()],
            MONTHS[usize::from(self.month.clamp(1, 12)) - 1],
            self.day,
            self.year
        )
    }

    /// A stamp from the original's two texts: a date "m/d/y" and a time "h:mm:ss".
    pub fn from_texts(date: &str, time: &str) -> Option<Stamp> {
        let numbers = |text: &str, sep: char| -> Option<Vec<u32>> {
            text.split(sep).map(|n| n.trim().parse().ok()).collect()
        };
        let [month, day, year] = <[u32; 3]>::try_from(numbers(date, '/')?).ok()?;
        let [hour, minute, second] = <[u32; 3]>::try_from(numbers(time, ':')?).ok()?;
        // two-digit years, as in "01/01/00", are in this century
        let year = if year < 100 { year + 2000 } else { year };
        Some(Stamp {
            year: u16::try_from(year).ok()?,
            month: u8::try_from(month).ok().filter(|m| (1..=12).contains(m))?,
            day: u8::try_from(day).ok().filter(|d| (1..=31).contains(d))?,
            hour: u8::try_from(hour).ok()?,
            minute: u8::try_from(minute).ok()?,
            second: u8::try_from(second).ok()?,
        })
    }
}

/// A Sanger trace: the signal in four channels and the base calls along it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Trace {
    pub sample: String,
    /// The instrument's model or name, as the file gives it.
    pub instrument: String,
    /// The base each channel stands for, in the order of `channels` (for example G, A, T, C).
    pub order: [u8; 4],
    /// The signal of each channel at every scan; the analysed signal when the file has it, else the raw one.
    pub channels: [Vec<i16>; 4],
    /// The bases the machine called.
    pub bases: Vec<u8>,
    /// Phred-style quality of each base; all zero from an instrument that gives none.
    pub quality: Vec<u8>,
    /// The scan at which each base's peak lies.
    pub peaks: Vec<u32>,
    /// How strong the signal of each channel is (`S/N%`), in the order of `channels`; zero when the file gives none.
    pub strength: [i16; 4],
    pub started: Option<Stamp>,
    pub ended: Option<Stamp>,
}

impl Trace {
    /// How many scans the traces have.
    pub fn scans(&self) -> usize {
        self.channels.iter().map(Vec::len).max().unwrap_or(0)
    }

    /// The highest signal anywhere in the trace.
    pub fn tallest(&self) -> i16 {
        self.channels
            .iter()
            .flatten()
            .copied()
            .max()
            .unwrap_or(0)
            .max(1)
    }

    /// Whether the machine gave quality values at all.
    pub fn has_quality(&self) -> bool {
        self.quality.iter().any(|&q| q > 0)
    }

    /// The share of bases with a quality of at least `q`, in percent.
    pub fn percent_at_least(&self, q: u8) -> f64 {
        if self.quality.is_empty() {
            return 0.0;
        }
        100.0 * self.quality.iter().filter(|&&x| x >= q).count() as f64 / self.quality.len() as f64
    }

    /// The channel (0 to 3) that stands for `base`.
    #[cfg(test)]
    pub fn channel_of(&self, base: u8) -> Option<usize> {
        self.order
            .iter()
            .position(|&b| b == base.to_ascii_uppercase())
    }

    /// The stretch `[start, end)` of bases worth keeping: the segment with the most reliable bases, found as
    /// Richard Mott's modified algorithm does it — each base scores `cutoff` minus its probability of being wrong
    /// (from its quality), and the stretch with the highest total wins. `None` when no stretch scores above zero
    /// (no quality values, or all poor).
    pub fn good_stretch(&self, cutoff: f64) -> Option<(usize, usize)> {
        let (mut best, mut best_range) = (0.0, None);
        let (mut running, mut from) = (0.0, 0usize);
        for (i, &q) in self.quality.iter().enumerate() {
            let score = cutoff - 10f64.powf(-f64::from(q) / 10.0);
            if running <= 0.0 {
                (running, from) = (0.0, i);
            }
            running += score;
            if running > best {
                best = running;
                best_range = Some((from, i + 1));
            }
        }
        best_range
    }
}

/// The run's start or end of number `number` (1 or 2): `RUND` holds the year (two bytes), month and day, `RUNT` the
/// hour, minute and second. A run that gives only one of them keeps the original's default for the other.
fn stamp(file: &Abif, number: i32) -> Option<Stamp> {
    let date = file.data(b"RUND", number).filter(|d| d.len() == 4);
    let time = file.data(b"RUNT", number).filter(|t| t.len() >= 3);
    if date.is_none() && time.is_none() {
        return None;
    }
    let mut s = Stamp::default();
    if let Some(d) = date {
        (s.year, s.month, s.day) = (u16::from_be_bytes([d[0], d[1]]), d[2], d[3]);
    }
    if let Some(t) = time {
        (s.hour, s.minute, s.second) = (t[0], t[1], t[2]);
    }
    Some(s)
}

/// Read a trace from the bytes of an `.ab1` file, by the rules of the original's reader (`read_abi_info`,
/// AppMain.tcl:30572): the base calls, peak positions and qualities of set 1 are the ones shown (the original calls
/// them "edited"; most files carry two identical sets), the signal is the analysed one (`DATA` 9 to 12), the channels
/// are in the order `FWO_` gives, or A, C, G, T when it gives none.
pub fn parse(bytes: &[u8]) -> Result<Trace, String> {
    let file = Abif::open(bytes)?;
    // set 2 is only a fallback, for a file that has no set 1
    let set = if file.data(b"PBAS", 1).is_some() {
        1
    } else {
        2
    };
    let bases = file
        .bytes(b"PBAS", set)
        .ok_or("The file has no base calls: it is not a sequencing trace.")?;
    let peaks = file.words(b"PLOC", set).unwrap_or_default();
    let mut quality = file.bytes(b"PCON", set).unwrap_or_default();
    let n = bases.len().min(if peaks.is_empty() {
        bases.len()
    } else {
        peaks.len()
    });
    quality.resize(n, 0);
    let mut order = *b"ACGT";
    if let Some(o) = file.text(b"FWO_", 1).filter(|o| o.len() == 4) {
        order.copy_from_slice(&o.to_ascii_uppercase().into_bytes());
    }
    let mut strength = [0i16; 4];
    if let Some(found) = file.shorts(b"S/N%", 1) {
        for (slot, value) in strength.iter_mut().zip(found) {
            *slot = value;
        }
    }
    let signal = |numbers: [i32; 4]| -> Option<[Vec<i16>; 4]> {
        let mut found = numbers.map(|k| file.shorts(b"DATA", k));
        found
            .iter()
            .all(Option::is_some)
            .then(|| found.each_mut().map(|c| c.take().unwrap_or_default()))
    };
    let channels = signal(ANALYSED)
        .or_else(|| signal(RAW))
        .ok_or("The file has no signal: the four colour traces are missing.")?;
    Ok(Trace {
        sample: file.text(b"SMPL", 1).unwrap_or_default(),
        instrument: file
            .text(b"MODL", 1)
            .filter(|m| !m.is_empty())
            .or_else(|| file.text(b"MCHN", 1))
            .unwrap_or_default(),
        order,
        channels,
        bases: bases[..n].to_vec(),
        quality,
        peaks: peaks.get(..n).map_or_else(
            || (0..n).map(|i| i as u32).collect(),
            |p| p.iter().map(|&x| u32::from(x)).collect(),
        ),
        strength,
        started: stamp(&file, 1),
        ended: stamp(&file, 2),
    })
}

/// Read a trace file of either kind, told by its first four bytes as the original does: ABIF (`.ab1`) or SCF
/// (`.scf`).
pub fn read(bytes: &[u8]) -> Result<Trace, String> {
    match bytes.get(..4) {
        Some(b"ABIF") => parse(bytes),
        Some(b".scf") => super::scf::parse(bytes),
        _ => Err("This is not a sequencing trace: it starts neither as an ABIF file (.ab1) nor as an SCF file (.scf).".into()),
    }
}

#[cfg(test)]
mod tests;
