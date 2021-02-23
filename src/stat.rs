//! Types to track information about number of matches and quality information
use std::cmp::Reverse;
use std::{collections::HashMap, io};

use ndarray::Array1;

use crate::match_type;
use crate::output::write_qual_report_header;
/// General information collected during read filtering
///
/// Usage: mutably borrowed by `read_filter::ReadFilter`
#[derive(Default, Debug)]
pub struct RunningStats {
    pub total_reads: u32,
    pub matching_reads: u32,
    pub ambigiuous_rejected: u32,
    pub peak_rejected: u32,
    pub mean_rejected: u32,
}

/// Tracker of important QC information
///
/// How are quality scores distributed accross the region of interest?
/// How accurate is the position information?
/// Append the records to this map.
pub struct QualStats {
    dat: HashMap<(u32, u8, u8), QualStatEntry>,
}

impl QualStats {
    pub fn new() -> Self {
        QualStats {
            dat: HashMap::new(),
        }
    }

    pub fn append(&mut self, mat: &match_type::SearchMatch) {
        let key = (mat.start_pos, mat.peak_qual(), mat.mean_qual());
        self.dat
            .entry(key)
            .and_modify(|existing| *existing += mat.into())
            .or_insert_with(|| mat.into());
    }
    fn entries_ordered(&self) -> QualStatsIter {
        let mut keys: Vec<_> = self.dat.keys().cloned().collect();
        keys.sort_unstable_by_key(|&e| Reverse(e)); // TODO: Verify behavior with tuple key
        QualStatsIter {
            ord_keys: keys,
            qs: self,
        }
    }
    #[allow(dead_code)]
    fn entries_unordered(&self) -> QualStatsIter {
        let keys: Vec<_> = self.dat.keys().into_iter().cloned().collect();
        QualStatsIter {
            ord_keys: keys,
            qs: self,
        }
    }

    /// Write the QC report directly as a .tsv
    ///
    /// By filtering by `dist_start`, `peak_qual`, `mean_qual` downstream tools can identify more appropriate filter values.
    pub fn write_to_buf<T: io::Write>(&self, buf: &mut T, seq_len: usize) -> io::Result<()> {
        write_qual_report_header(buf, seq_len)?;
        for (k, v) in self.entries_ordered() {
            write!(
                buf,
                "{dist_start}\t{peak_qual}\t{mean_qual}\t{reads}\t{reverse_reads}",
                dist_start = k.0,
                peak_qual = k.1,
                mean_qual = k.2,
                reads = v.reads(),
                reverse_reads = v.reverse_reads()
            )?;
            for e in v.normalized_qual().iter() {
                write!(buf, "\t{}", e)?;
            }
            write!(buf, "\n")?;
        }
        Ok(())
    }
}

struct QualStatsIter<'a> {
    ord_keys: Vec<(u32, u8, u8)>,
    qs: &'a QualStats,
}

impl<'a> Iterator for QualStatsIter<'a> {
    type Item = ((u32, u8, u8), &'a QualStatEntry);

    fn next(&mut self) -> Option<Self::Item> {
        let k = self.ord_keys.pop()?;
        let v = self.qs.dat.get(&k)?;
        Some((k, v))
    }
}

struct QualStatEntry {
    read_count: u32,
    reverse_count: u32,
    qual_arr: Array1<u32>,
}

impl QualStatEntry {
    fn normalized_qual(&self) -> Array1<f32> {
        let arr: Array1<f32> = self.qual_arr.mapv(|i| i as f32);
        arr / (self.read_count as f32) - 33.0
    }
    fn reads(&self) -> u32 {
        self.read_count
    }
    fn reverse_reads(&self) -> u32 {
        self.reverse_count
    }
}

impl From<&match_type::SearchMatch> for QualStatEntry {
    // TODO: Decide if this is smart given the always known seq_len
    fn from(mat: &match_type::SearchMatch) -> Self {
        let qual_arr: Array1<u32> = mat.quality.iter().map(|&a| a as u32).collect();
        QualStatEntry {
            read_count: 1,
            reverse_count: mat.reverse_strand as u32,
            qual_arr,
        }
    }
}

impl std::ops::AddAssign for QualStatEntry {
    fn add_assign(&mut self, rhs: Self) {
        self.read_count += rhs.read_count;
        self.reverse_count += rhs.reverse_count;
        self.qual_arr += &rhs.qual_arr; // Would panic if shapes mismatch
    }
}
