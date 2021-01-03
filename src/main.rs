use bio::alphabets::dna;
use bio::io::fastq;
use counter::Counter;
use fastq::Records;
use read_filter::handling::{GracefulOption, GracefulResult};
use read_filter::jsonconf;
use regex::bytes::Regex;
use std::fmt;
use std::iter::Iterator;
#[allow(unused_imports)]
use std::{error::Error, fmt::Display, todo};

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
    // Same as before...
    let infile = std::path::Path::new(&cfg.infile);
    let bname = infile
        .file_name()
        .unwrap_graceful("Input needs to be a file"); // TODO file name is just basename only errs when given ..
    let oname = String::from(bname.to_str().unwrap())
        .replace(".fastq.gz", OUT_ENDING)
        .replace(".txt.gz", OUT_ENDING)
        .replace(".fastq", OUT_ENDING);
    let outdir = std::path::Path::new(&cfg.outdir);
    let outfile = outdir.join(oname);

    // FASTQ parsing
    let (reader, _) = niffler::from_path(infile).unwrap_formatful("Invalid input path!");
    let fq_reader = fastq::Reader::new(reader);

    let mut stats = RunningStats::default();
    let rf = ReadFilter::new(fq_reader.records(), &cfg, &mut stats);

    let c = rf
        .map(move |a| String::from_utf8(a.content).unwrap())
        .collect::<Counter<_>>();
    println!("{:?}", stats)
    // println!("{:?}", c);
    // for sm in rf {
    //     println!("{}, Peak {}, Mean {}", sm, sm.peak_qual(), sm.mean_qual());
    // }
}

#[derive(Debug, PartialEq)]
struct SearchMatch {
    // TODO: lazy alloc alternative that only allocs if a match is confirmed to be kept
    // Currently we don't care about types for optimal memory foot print or alignment. Simpler ownership model is preferred over optimization of heap allocations.
    // Array of structs might be more practical than struct of arrays as the primary application is streaming intensive
    content: Vec<u8>,
    // We need to own the content at some point
    quality: Vec<u8>,
    // Keeping the quality string around may be optional.
    // inherent properties like mean or peak qual score can be computed trivially or alternatively stored if memory footprind would be a serious concern.
    reverse_strand: bool,
    // Relevant for reporting
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
    let regex = format!("(?-u)^.{{{expt_begin},{expt_end}}}{left_flank}([ACGT]{{{content_length}}}){right_flank}.*$",
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
struct ReadFilter<'a, T: std::io::Read> {
    fq_records: Records<T>,
    regex: Regex,
    min_mean_qual: Option<u8>,
    min_peak_qual: Option<u8>,
    stats: &'a mut RunningStats,
}

impl<'a, T> ReadFilter<'a, T>
where
    T: std::io::Read,
{
    fn new(fq_records: Records<T>, cfg: &ProgConfig, stats: &'a mut RunningStats) -> Self {
        let expt_begin = cfg.expected_start.saturating_sub(cfg.position_tolerance);
        let expt_end = cfg.expected_start + cfg.position_tolerance;
        let regex = format!("(?-u)^.{{{expt_begin},{expt_end}}}{left_flank}([ACGT]{{{content_length}}}){right_flank}.*$",
            expt_begin=expt_begin,
            expt_end=expt_end,
            left_flank=cfg.left_flank,
            right_flank=cfg.right_flank,
            content_length=cfg.insert_length);
        let regex = Regex::new(&regex).expect("Error while compiling regex"); // Any sensible checks that should be done before?

        let min_mean_qual = cfg.min_mean_qual;
        let min_peak_qual = cfg.min_peak_qual;

        ReadFilter {
            fq_records,
            regex,
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
            let rec = self.fq_records.next()?.ok()?;
            self.stats.total_reads += 1;
            let provisional = match match_both_strands(&rec, &self.regex) {
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

                break result;
            }
        };
        Some(res)
    }
}

// TODO Try out regex free alternative
fn match_both_strands(
    read: &bio::io::fastq::Record,
    compiled_expression: &Regex,
) -> (Option<SearchMatch>, Option<SearchMatch>) {
    let rev_seq = dna::revcomp(read.seq()); // TODO: Remove eager alloc
    // TODO: Reduce Code duplication
    let fwd = if let Some(captures) = compiled_expression.captures(read.seq()) {
        let content_match = captures
            .get(1)
            .expect("capture group should contain the sequence content");
        let content = content_match.as_bytes().to_vec();
        let start_pos = content_match.start();
        let quality: Vec<u8> = (&read.qual()[start_pos..content_match.end()]).to_vec();
        Some(SearchMatch {
            content,
            quality,
            reverse_strand: false,
            start_pos: (start_pos as u32),
        })
    } else {
        None
    };

    let rev = if let Some(captures) = compiled_expression.captures(&rev_seq) {
        let content_match = captures
            .get(1)
            .expect("capture group should contain the sequence content");
        let content = content_match.as_bytes().to_vec();
        let start_pos = content_match.start();
        // Reversing quality
        let end_pos = content_match.end();
        let seq_len = read.seq().len();
        let mut quality: Vec<u8> = (&read.qual()[seq_len - end_pos..seq_len - start_pos]).to_vec();
        quality.reverse();
        Some(SearchMatch {
            content,
            quality,
            reverse_strand: true,
            start_pos: (start_pos as u32),
        })
    } else {
        None
    };
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
            min_peak_qual: json_config.qual_peak, // TODO: Maybe change the config type to accept options
            min_mean_qual: json_config.qual_mean,
        })
    }
}
