use crate::config::ProgConfig;
use crate::match_type::{CandidateMatch, SearchMatch};
use crate::stat::RunningStats;
use crate::ExactPattern;
use bio::alphabets::dna;
use bio::io::fastq;
use fastq::Records;
use std::io;

pub fn is_dna_char(char: &u8) -> bool {
    matches!(char, b'A' | b'C' | b'G' | b'T')
}

pub struct PrecomputedPatterns {
    pub fwd_start: ExactPattern,
    pub fwd_end: Vec<u8>,
    pub rev_start: Vec<u8>,
    pub rev_end: ExactPattern,
    pub content_len: usize,
    pub start_len: usize,
    pub end_len: usize,
    pub expt_begin: usize,
    pub expt_end: usize,

    pub fwd_dist: usize,
    pub rev_dist: usize,
    pub total_len: usize,
}

impl PrecomputedPatterns {
    pub fn new(cfg: &ProgConfig, expt_begin: usize, expt_end: usize) -> Self {
        let start_len = cfg.left_flank.len();
        let end_len = cfg.right_flank.len();
        let content_len = cfg.insert_length as usize;
        let total_len = start_len + content_len + end_len;
        let fwd_dist = start_len + content_len;
        let rev_dist = end_len + content_len;
        // Performs the String to Vec<u8> cast
        let fwd_start = ExactPattern::new(cfg.left_flank.as_bytes());
        let fwd_end = cfg.right_flank.as_bytes().to_vec();
        let rev_start = dna::revcomp(cfg.left_flank.as_bytes());
        let rev_end = ExactPattern::new(dna::revcomp(cfg.right_flank.as_bytes()));

        PrecomputedPatterns {
            fwd_start,
            fwd_end,
            rev_start,
            rev_end,
            content_len,
            start_len,
            end_len,
            expt_begin,
            expt_end,
            fwd_dist,
            rev_dist,
            total_len,
        }
    }
}

pub struct ReadFilter<'a, T: io::Read> {
    fq_records: Records<T>,
    pats: PrecomputedPatterns,
    min_mean_qual: Option<u8>,
    min_peak_qual: Option<u8>,
    stats: &'a mut RunningStats,
}

impl<'a, T> ReadFilter<'a, T>
where
    T: io::Read,
{
    pub fn new(fq_records: Records<T>, cfg: &ProgConfig, stats: &'a mut RunningStats) -> Self {
        let expt_begin = cfg.expected_start.saturating_sub(cfg.position_tolerance) as usize;
        let expt_end = (cfg.expected_start + cfg.position_tolerance) as usize;

        let min_mean_qual = cfg.min_mean_qual;
        let min_peak_qual = cfg.min_peak_qual;
        let pats = PrecomputedPatterns::new(cfg, expt_begin, expt_end);

        ReadFilter {
            fq_records,
            pats,
            min_mean_qual,
            min_peak_qual,
            stats,
        }
    }
}

impl<'a, T> Iterator for ReadFilter<'a, T>
where
    T: io::BufRead,
{
    type Item = SearchMatch;

    fn next(&mut self) -> Option<Self::Item> {
        let res = loop {
            let rec = self.fq_records.next()?.ok()?; // ParseErrors are silenced!
            self.stats.total_reads += 1;
            let provisional = match match_both_strands(&rec, &self.pats) {
                (Some(a), None) => Some(a),
                (None, Some(b)) => Some(b),
                (Some(_), Some(_)) => {
                    self.stats.ambigiuous_rejected += 1;
                    None
                } // TODO: How to deal with this illegal situation?
                _ => None,
            };
            if let Some(result) = provisional {
                if !result.seq().iter().all(is_dna_char) {
                    continue;
                }
                self.stats.matching_reads += 1;
                if let Some(min) = self.min_peak_qual {
                    if result.peak_qual() < min {
                        self.stats.peak_rejected += 1;
                        continue;
                    }
                }
                if let Some(min) = self.min_mean_qual {
                    if result.mean_qual() < min {
                        self.stats.mean_rejected += 1;
                        continue;
                    }
                }

                // First heap allocs after fastq parse happen here
                break result.materialize();
            } // Else loop again till match or exhaustion
        };
        Some(res)
    }
}

pub fn match_both_strands<'a>(
    read: &'a bio::io::fastq::Record,
    patterns: &PrecomputedPatterns,
) -> (Option<CandidateMatch<'a>>, Option<CandidateMatch<'a>>) {
    let read_seq = read.seq();
    let read_len = read_seq.len();

    let mat_fwd = patterns.fwd_start.find_all(read_seq).find(|&idx| {
        (idx>= patterns.expt_begin)
         && (idx<= patterns.expt_end)
         && (idx+patterns.total_len<= read_len) // Necessary to ensure legal indexing
         && (read_seq[idx+patterns.fwd_dist..idx+patterns.total_len]==patterns.fwd_end)
    });

    let fwd = if let Some(idx) = mat_fwd {
        let start_idx = idx + patterns.start_len;
        let range = start_idx..start_idx + patterns.content_len;
        let mat = CandidateMatch::new(
            &read_seq[range.clone()],
            &read.qual()[range],
            false,
            start_idx as u32,
        );
        Some(mat)
    } else {
        None
    };

    let mat_rev = patterns.rev_end.find_all(read_seq).find(
        |&idx| {
            (idx + patterns.total_len + patterns.expt_begin <= read_len)
                && (idx + patterns.total_len + patterns.expt_end >= read_len)
                && (read_seq[idx + patterns.rev_dist..idx + patterns.total_len]
                    == patterns.rev_start)
        }, // Legal due to first condition
    );
    let rev = if let Some(idx) = mat_rev {
        let start_pos = read_len - (idx + patterns.rev_dist); // No underflow possible due to first condition
        let range = idx + patterns.end_len..idx + patterns.rev_dist;
        let mat = CandidateMatch::new(
            &read_seq[range.clone()],
            &read.qual()[range],
            true,
            start_pos as u32,
        );
        Some(mat)
    } else {
        None
    };

    (fwd, rev)
}
