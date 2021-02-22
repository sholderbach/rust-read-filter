//! Helpers to write output to buffers

use std::io;

use crate::config::ProgConfig;
use crate::stat::RunningStats;

pub fn write_stats_header<T: io::Write>(buf: &mut T, stats: &RunningStats) -> io::Result<()> {
    let qual_reads = stats.matching_reads - (stats.peak_rejected + stats.mean_rejected);
    write!(
        buf,
        "# raw_total_reads: {total_reads}\n\
            # matching_reads: {matching_reads}\n\
            # peak_qual_rejected_reads: {peak_rejected}\n\
            # mean_qual_rejected_reads: {mean_rejected}\n\
            # ambiguous_matches_rejected: {ambiguous_rejected}\n\
            # quality_reads: {qual_reads}\n",
        total_reads = stats.total_reads,
        matching_reads = stats.matching_reads,
        peak_rejected = stats.peak_rejected,
        mean_rejected = stats.mean_rejected,
        ambiguous_rejected = stats.ambigiuous_rejected,
        qual_reads = qual_reads,
    )
}

pub fn write_config_header<T: io::Write>(buf: &mut T, cfg: &ProgConfig) -> io::Result<()> {
    // Writing the regex is to reflect the original python version, but no guarantee that we use the exact regex
    // Code duplication or out of sync risk
    let expt_begin = cfg.expected_start.saturating_sub(cfg.position_tolerance);
    let expt_end = cfg.expected_start + cfg.position_tolerance;
    let regex = format!(
        "^.{{{expt_begin},{expt_end}}}{left_flank}([ACGT]{{{content_length}}}){right_flank}.*$",
        expt_begin = expt_begin,
        expt_end = expt_end,
        left_flank = cfg.left_flank,
        right_flank = cfg.right_flank,
        content_length = cfg.insert_length
    );
    write!(
        buf,
        "# filter: {regex}\n\
            # accepted_peak_qual: {peak_accepted}\n\
            # accepted_mean_qual: {mean_accepted}\n",
        regex = regex,
        peak_accepted = cfg.min_peak_qual.unwrap_or_default(),
        mean_accepted = cfg.min_mean_qual.unwrap_or_default(),
    )
}

pub fn write_read_report_header<T: io::Write>(buf: &mut T, seq_len: usize) -> io::Result<()> {
    write!(buf, "read\tdist_start\treversed\tpeak_qual\tmean_qual")?;
    for i in 0..seq_len {
        write!(buf, "\tqual_pos_{}", i)?;
    }
    write!(buf, "\n")
}

pub fn write_qual_report_header<T: io::Write>(buf: &mut T, seq_len: usize) -> io::Result<()> {
    write!(
        buf,
        "dist_start\tpeak_qual\tmean_qual\treads\treverse_reads"
    )?;
    for i in 0..seq_len {
        write!(buf, "\tqual_pos_{}", i)?;
    }
    write!(buf, "\n")
}
