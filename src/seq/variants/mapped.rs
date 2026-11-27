//! Mapping reads to a reference with the program's mapper (`seq::minimap`, with the base-level alignment) and keeping the mappings as SAM
//! records: every read gives its alignments, or one line that says it did not map. A read may lie in several places that do not overlap on it
//! (a chimera, or a read that crosses the start of a circle): SAM has one primary line for a read, so the best is the primary and each of the
//! others is flagged supplementary (2048). `NH` is the number of lines for one place: the primary or supplementary line and its secondaries.

use super::reads::Read;
use super::sam::{Record, Sam, SUPPLEMENTARY};
use crate::seq::minimap::{map_read, Index, Named, Options};
use crate::seq::par::par_map;

/// The mappings of `reads` to `targets`, in the order of the reads. `Err` when no target can be indexed, or when `o` does not ask for the
/// alignment of the bases (`cigar`): a SAM line has the CIGAR.
pub fn map_to_sam(targets: &[Named], reads: &[Read], o: &Options) -> Result<Sam, String> {
    if !o.cigar {
        return Err(
            "SAM needs the alignment of the bases: map with the base-level alignment on (-c)."
                .into(),
        );
    }
    if o.overlaps {
        return Err("Overlaps between reads are PAF, not SAM: a SAM line lies on a target.".into());
    }
    let index = Index::build(targets.to_vec(), o)?;
    let ids: Vec<usize> = (0..reads.len()).collect();
    let rows = par_map(&ids, |&i| {
        let read = &reads[i];
        let hits = map_read(&index, &read.seq, o, &|_| false);
        let mut records: Vec<Record> = Vec::new();
        // the mapper gives each place as its primary followed by its secondaries, the best place first
        for place in hits.chunk_by(|_, next| !next.primary) {
            for h in place {
                if let Some(mut record) = Record::from_hit(read, &index, h) {
                    if h.primary && !records.is_empty() {
                        record.flag |= SUPPLEMENTARY;
                    }
                    record.tags.push(format!("NH:i:{}", place.len()));
                    records.push(record);
                }
            }
        }
        if records.is_empty() {
            vec![Record::unmapped(read)]
        } else {
            records
        }
    });
    Ok(Sam {
        targets: targets
            .iter()
            .map(|(name, bases)| (name.clone(), bases.len()))
            .collect(),
        records: rows.into_iter().flatten().collect(),
    })
}

#[cfg(test)]
mod tests;
