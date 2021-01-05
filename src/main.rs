use bio::alphabets::dna;
use bio::io::fastq;
use counter::Counter;
use fastq::Records;
use read_filter::handling::{GracefulOption, GracefulResult};
use read_filter::jsonconf;
use std::{fmt, io::BufWriter, io::Write};
use std::iter::Iterator;
#[allow(unused_imports)]
use std::{error::Error, fmt::Display, todo};
#[allow(unused_imports)]
use bio::pattern_matching::{shift_and, bndm};

type ExactPattern = shift_and::ShiftAnd;

#[macro_use]
extern crate clap;

const OUT_ENDING: &str = ".processed.tsv";
#[allow(unused)]
const RR_ENDING: &str = ".readreport.tsv";
#[allow(unused)]
const QC_ENDING: &str = ".quality.tsv";
#[allow(unused)]
fn main() {
    let cfg = ProgConfig::from_cli().unwrap_graceful();
    let infile = std::path::Path::new(&cfg.infile);
    let bname = infile
        .file_name()
        .unwrap_graceful("Input needs to be a file"); // TODO file name is just basename only errs when given ..
    let oname = String::from(bname.to_str().unwrap())
        .replace(".fastq.gz", OUT_ENDING)
        .replace(".txt.gz", OUT_ENDING)
        .replace(".fastq", OUT_ENDING);
    let outdir = std::path::Path::new(&cfg.outdir);
    std::fs::create_dir_all(outdir).unwrap_messageful(&format!("Could not create output directory at: {:?}", outdir.to_str().unwrap()));
    let outfile = outdir.join(oname);
    let ofile = std::fs::File::create(&outfile).unwrap_messageful(&format!("Could not create output file at: {:?}", outfile.to_str().unwrap()));

    // FASTQ parsing
    let (reader, _) = niffler::from_path(infile).unwrap_formatful("Invalid input path!");
    let fq_reader = fastq::Reader::new(reader);

    let mut stats = RunningStats::default();
    let rf = ReadFilter::new(fq_reader.records(), &cfg, &mut stats);

    let c = rf
        .map(move |a| a.content)
        .collect::<Counter<_>>();
    
    let mut ofile= BufWriter::new(ofile);
    write_config_header(&mut ofile, &cfg).unwrap_messageful("Error while writing output");
    write_stats_header(&mut ofile, &stats);
    writeln!(ofile, "seq\treads");
    for (seq, count) in c.iter() {
        writeln!(ofile, "{}\t{}", std::str::from_utf8(seq).unwrap(), count);
    }
}

#[derive(Debug, PartialEq)]
struct SearchMatch {
    content: Vec<u8>,
    quality: Vec<u8>, // TODO maybe have a variant that doesn't keep quality if not needed
    reverse_strand: bool,
    start_pos: u32,
    // POTENTIAL FUTURE
    // pattern_id: u32,
}

impl SearchMatch {
    fn peak_qual(&self) -> u8 {
        self.quality
            .iter()
            .min()
            .expect("Expect the extraction of nonempty content")
            - 33
    }
    fn mean_qual(&self) -> u8 {
        let avg_qual =
            self.quality.iter().fold(0u32, |x, b| x + (*b as u32)) / self.quality.len() as u32 - 33;
        avg_qual as u8
    }

    fn to_fastq(&self, id: &str) -> fastq::Record {
        let strand = if self.reverse_strand { "-" } else { "+" };
        fastq::Record::with_attrs(id, Some(strand), &self.content, &self.quality)
    }

    fn to_fastq_original_strand(&self, id: &str) -> fastq::Record {
        let strand = if self.reverse_strand { "-" } else { "+" };
        if !self.reverse_strand {
            fastq::Record::with_attrs(id, Some(strand), &self.content, &self.quality)
        } else {
            fastq::Record::with_attrs(
                id,
                Some(strand),
                &dna::revcomp(&self.content),
                &self.quality.iter().rev().copied().collect::<Vec<_>>(),
            )
        }
    }
}
impl Display for SearchMatch {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let strand = if self.reverse_strand { "-" } else { "+" };
        let seq: &str = std::str::from_utf8(&self.content).unwrap_or_default();
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

#[derive(Default, Debug)]
struct RunningStats {
    total_reads: u32,
    matching_reads: u32,
    ambigiuous_rejected: u32,
    peak_rejected: u32,
    mean_rejected: u32,
}

fn write_stats_header<T: std::io::Write>(buf: &mut T, stats: &RunningStats) -> std::io::Result<()> {
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

fn write_config_header<T: std::io::Write>(buf: &mut T, cfg: &ProgConfig) -> std::io::Result<()> {
    // Writing the regex is to reflect the original python version, but no guarantee that we use the exact regex
    // Code duplication or out of sync risk
    let expt_begin = cfg.expected_start.saturating_sub(cfg.position_tolerance);
    let expt_end = cfg.expected_start + cfg.position_tolerance;
    let regex = format!("^.{{{expt_begin},{expt_end}}}{left_flank}([ACGT]{{{content_length}}}){right_flank}.*$",
            expt_begin=expt_begin,
            expt_end=expt_end,
            left_flank=cfg.left_flank,
            right_flank=cfg.right_flank,
            content_length=cfg.insert_length);
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


struct PrecomputedPatterns {
    fwd_start: ExactPattern,
    fwd_end: Vec<u8>,
    rev_start: Vec<u8>, 
    rev_end: ExactPattern,
    content_len: usize,
    start_len: usize,
    end_len:usize,
    expt_begin: usize,
    expt_end: usize,

    fwd_dist: usize,
    rev_dist: usize,
    total_len: usize,
}

struct ReadFilter<'a, T: std::io::Read> {
    fq_records: Records<T>,
    pats: PrecomputedPatterns,
    min_mean_qual: Option<u8>,
    min_peak_qual: Option<u8>,
    stats: &'a mut RunningStats,
}

impl<'a, T> ReadFilter<'a, T>
where
    T: std::io::Read,
{
    fn new(fq_records: Records<T>, cfg: &ProgConfig, stats: &'a mut RunningStats) -> Self {
        let expt_begin = cfg.expected_start.saturating_sub(cfg.position_tolerance) as usize;
        let expt_end = (cfg.expected_start + cfg.position_tolerance) as usize;

        let min_mean_qual = cfg.min_mean_qual;
        let min_peak_qual = cfg.min_peak_qual;

        let start_len = cfg.left_flank.len();
        let end_len = cfg.right_flank.len();
        let content_len= cfg.insert_length as usize;
        let total_len = start_len + content_len + end_len;
        let fwd_dist = start_len + content_len;
        let rev_dist = end_len + content_len;
        let fwd_start = ExactPattern::new(cfg.left_flank.as_bytes());
        let fwd_end = cfg.right_flank.as_bytes().to_vec();
        let rev_start = dna::revcomp(cfg.left_flank.as_bytes());
        let rev_end = ExactPattern::new(dna::revcomp(cfg.right_flank.as_bytes()));

        let pats = PrecomputedPatterns{
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
        };

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
    T: std::io::Read,
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

/// Zero-copy version of a match for internal processing
/// the raw `seq` and `qual` fields are not adjusted to account for the correct strand
struct CandidateMatch<'a> {
    seq: &'a [u8],
    quality: &'a [u8],
    reverse_strand: bool,
    start_pos: u32, // Adjusted to common strandedness? TODO
}

impl<'a> CandidateMatch<'a>{
    fn materialize(self) -> SearchMatch {
        //! Consumes self to produce an owned `SearchMatch`
        //! calls `dna::revcomp` to produce the useful sequence orientation
        if !self.reverse_strand {
            SearchMatch {
                content: self.seq.to_vec(),
                quality: self.quality.to_vec(),
                reverse_strand: self.reverse_strand,
                start_pos: self.start_pos,
            }
        } else {
            let mut quality = self.quality.to_vec();
            quality.reverse();
            SearchMatch {
                content: dna::revcomp(self.seq),
                quality,
                reverse_strand: self.reverse_strand,
                start_pos: self.start_pos,
            }
        }
    }

    fn peak_qual(&self) -> u8 {
        //! Lowest quality score found
        self.quality
            .iter()
            .min()
            .expect("Expect the extraction of nonempty content")
            - 33
    }
    fn mean_qual(&self) -> u8 {
        //! Average quality score
        //! Performs integer division
        let avg_qual =
            self.quality.iter().fold(0u32, |x, b| x + (*b as u32)) / self.quality.len() as u32 - 33;
        // SAFETY: as u8 considered safe as valid PHRED string assumed
        avg_qual as u8
    }
}

fn match_both_strands<'a>(
    read: &'a bio::io::fastq::Record,
    patterns: &PrecomputedPatterns,
) -> (Option<CandidateMatch<'a>>, Option<CandidateMatch<'a>>) {
    let read_seq = read.seq();
    let read_len = read_seq.len();

    let mat_fwd = patterns.fwd_start.find_all(read_seq).filter(|&idx| 
        (idx>= patterns.expt_begin) 
         && (idx<= patterns.expt_end) 
         && (idx+patterns.total_len<= read_len) // Necessary to ensure legal indexing
         && (&read_seq[idx+patterns.fwd_dist..idx+patterns.total_len]==patterns.fwd_end)
        ).next();
    
    let fwd = if let Some(idx) = mat_fwd {
        let start_idx = idx + patterns.start_len;
        let range =start_idx..start_idx+patterns.content_len;
        let mat = CandidateMatch{
            seq: &read_seq[range.clone()],
            quality: &read.qual()[range],
            reverse_strand: false,
            start_pos: start_idx as u32,
        };
        Some(mat)
    } else {None};
    

    let mat_rev = patterns.rev_end.find_all(read_seq).filter(|&idx| 
        (idx + patterns.total_len + patterns.expt_begin <= read_len) 
         && (idx + patterns.total_len + patterns.expt_end >= read_len) 
         && (&read_seq[idx+patterns.rev_dist..idx+patterns.total_len]==patterns.rev_start) // Legal due to first condition
        ).next();
    let rev = if let Some(idx) = mat_rev {
        let start_pos = read_len - (idx + patterns.rev_dist); // No underflow possible due to first condition
        let range = idx+patterns.end_len..idx + patterns.rev_dist;
        let mat = CandidateMatch{
            seq: &read_seq[range.clone()],
            quality: &read.qual()[range],
            reverse_strand: true,
            start_pos: start_pos as u32,
        };
        Some(mat)
    } else {None};
    
    (fwd, rev)
}

#[allow(unused)]
struct ProgConfig {
    infile: String,
    outdir: String,
    rr_required: bool,
    qc_required: bool,
    left_flank: String,
    right_flank: String,
    insert_length: u32,
    expected_start: u32,
    position_tolerance: u32,
    min_peak_qual: Option<u8>,
    min_mean_qual: Option<u8>,
}

impl ProgConfig {
    fn from_cli() -> Result<ProgConfig, Box<dyn Error>> {
        // Specify CLI
        let matches = clap_app!(readfilter =>
        (version: "0.1")
        (author: "Stefan Holderbach")
        (about: "Read filter for amplicon sequencing with a defined region")
        (@arg CONFIG: -c --config +takes_value +required "Sets a custom config file")
        (@arg INPUT: +required "Sets the input file to use")
        (@arg OUTPUT: +required "Sets the output path")
        (@arg read_report: -r --("read-report") "Also output a table with QC information for each read")
        (@arg qc_report: -q --("qc-report") "Also output a table with overall QC information")
        (@arg debug: -d ... "Sets the level of debugging information")
    ).get_matches();
        // Unpack arguments
        let infile = matches
            .value_of("INPUT")
            .unwrap_graceful("Missing inputfile");
        let outdir = matches
            .value_of("OUTPUT")
            .unwrap_graceful("Missing output directory");
        let config_file = matches
            .value_of("CONFIG")
            .unwrap_graceful("Missing config file");
        let rr_required = matches.is_present("read_report");
        let qc_required = matches.is_present("qc_report");

        let json_config =
            jsonconf::load_json_config(config_file).unwrap_formatful("While parsing configuration");

        Ok(ProgConfig {
            infile: infile.to_string(),
            outdir: outdir.to_string(),
            rr_required,
            qc_required,
            left_flank: json_config.left_flank,
            right_flank: json_config.right_flank,
            insert_length: json_config.content_length,
            expected_start: json_config.expect_begin,
            position_tolerance: json_config.tolerance,
            min_peak_qual: json_config.qual_peak,
            min_mean_qual: json_config.qual_mean,
        })
    }
}
