//! Types to track information about number of matches and quality information

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
