//! Types holding an individual match
use bio::alphabets::dna;
use bio::io::fastq;
use std::fmt;
use std::fmt::Display;
use std::io;

/// Owned version of a succesful match
#[derive(Debug, PartialEq)]
pub struct SearchMatch {
    /// Sequence as bytes
    pub seq: Vec<u8>,
    /// Quality string
    ///
    /// ATTENTION:
    /// Direct encoding from fastq. For PHRED scores subtract 33!
    pub quality: Vec<u8>, // TODO maybe have a variant that doesn't keep quality if not needed
    /// Whether match occurred on the reverse complement strand
    pub reverse_strand: bool,
    /// Index from the (reverse complement) start starting the content sequence
    pub start_pos: u32,
}

impl SearchMatch {
    /// Lowest quality observed in the content of the match. Corresponds to a peak in error probability
    pub fn peak_qual(&self) -> u8 {
        self.quality
            .iter()
            .min()
            .expect("Expect the extraction of nonempty content")
            - 33
    }
    /// integer average quality in the content of the match. Equivalent to `CandidateMatch.mean_qual` used to filter reads.
    pub fn mean_qual(&self) -> u8 {
        let avg_qual =
            self.quality.iter().fold(0u32, |x, b| x + (*b as u32)) / self.quality.len() as u32 - 33;
        avg_qual as u8
    }
    /// Floating point average quality (Higher accuracy for diagnostics)
    pub fn accurate_mean_qual(&self) -> f32 {
        let avg_qual = self.quality.iter().fold(0u32, |x, b| x + (*b as u32)) as f32
            / self.quality.len() as f32
            - 33f32;
        avg_qual
    }

    /// By providing a sequence `id` a match can be turned into a `fastq::Record` with the sequence in the orientation defined by the pattern
    pub fn to_fastq(&self, id: &str) -> fastq::Record {
        let strand = if self.reverse_strand { "-" } else { "+" };
        fastq::Record::with_attrs(id, Some(strand), &self.seq, &self.quality)
    }

    /// By providing a sequence `id` a match can be turned into a `fastq::Record` with the sequence in the original orientation
    pub fn to_fastq_original_strand(&self, id: &str) -> fastq::Record {
        let strand = if self.reverse_strand { "-" } else { "+" };
        if !self.reverse_strand {
            fastq::Record::with_attrs(id, Some(strand), &self.seq, &self.quality)
        } else {
            fastq::Record::with_attrs(
                id,
                Some(strand),
                &dna::revcomp(&self.seq),
                &self.quality.iter().rev().copied().collect::<Vec<_>>(),
            )
        }
    }

    /// Output a single tab-separated record for diagnostics
    pub fn write_read_report_line<T: io::Write>(&self, buf: &mut T) -> io::Result<()> {
        let seq: &str = std::str::from_utf8(&self.seq).unwrap_or_default();
        write!(
            buf,
            "{seq}\t{dist_start}\t{rev}\t{peak_qual}\t{mean_qual}",
            seq = seq,
            dist_start = self.start_pos,
            rev = self.reverse_strand,
            peak_qual = self.peak_qual(),
            mean_qual = self.accurate_mean_qual(),
        )?;
        for &val in self.quality.iter() {
            write!(buf, "\t{}", val - 33)?;
        }
        write!(buf, "\n")
    }
}
impl Display for SearchMatch {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let strand = if self.reverse_strand { "-" } else { "+" };
        let seq: &str = std::str::from_utf8(&self.seq).unwrap_or_default();
        let qual: &str = std::str::from_utf8(&self.quality).unwrap_or_default();
        write!(
            f,
            "{pos}({strd}): {seq}, PHRED: {qual}",
            pos = self.start_pos,
            strd = strand,
            seq = seq,
            qual = qual
        )
    }
}

/// Zero-copy version of a match for internal processing
/// the raw `seq` and `qual` fields are not adjusted to account for the correct strand
pub struct CandidateMatch<'a> {
    seq: &'a [u8],
    quality: &'a [u8],
    reverse_strand: bool,
    start_pos: u32, // Adjusted to common strandedness? TODO
}

impl<'a> CandidateMatch<'a> {
    pub fn new(seq: &'a [u8], quality: &'a [u8], reverse_strand: bool, start_pos: u32) -> Self {
        CandidateMatch {
            seq,
            quality,
            reverse_strand,
            start_pos,
        }
    }

    pub fn materialize(self) -> SearchMatch {
        //! Consumes self to produce an owned `SearchMatch`
        //! calls `dna::revcomp` to produce the useful sequence orientation
        if !self.reverse_strand {
            SearchMatch {
                seq: self.seq.to_vec(),
                quality: self.quality.to_vec(),
                reverse_strand: self.reverse_strand,
                start_pos: self.start_pos,
            }
        } else {
            let mut quality = self.quality.to_vec();
            quality.reverse();
            SearchMatch {
                seq: dna::revcomp(self.seq),
                quality,
                reverse_strand: self.reverse_strand,
                start_pos: self.start_pos,
            }
        }
    }

    pub fn peak_qual(&self) -> u8 {
        //! Lowest quality score found
        self.quality
            .iter()
            .min()
            .expect("Expect the extraction of nonempty content")
            - 33
    }
    pub fn mean_qual(&self) -> u8 {
        //! Average quality score
        //! Performs integer division
        let avg_qual =
            self.quality.iter().fold(0u32, |x, b| x + (*b as u32)) / self.quality.len() as u32 - 33;
        // SAFETY: as u8 considered safe as valid PHRED string assumed
        avg_qual as u8
    }

    pub fn seq(&self) -> &'a [u8] {self.seq}
}
