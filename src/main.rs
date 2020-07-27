use bio::io::fastq;
use read_filter::handling::{GracefulOption, GracefulResult};
use read_filter::jsonconf;
use std::error::Error;
// TODO: Decide on json crate for config file

#[macro_use]
extern crate clap;

#[allow(unused)]
fn main() {
    const OUT_ENDING: &str = ".processed.tsv";
    const RR_ENDING: &str = ".readreport.tsv";
    const QC_ENDING: &str = ".quality.tsv";

    let cfg = ProgConfig::from_cli().unwrap_graceful();
    // Same as before...
    let infile = std::path::Path::new(&cfg.infile);
    let bname = infile
        .file_name()
        .unwrap_graceful("Input needs to be a file");
    let oname = String::from(bname.to_str().unwrap()) // Can this be considered a given?
        .replace(".fastq.gz", OUT_ENDING)
        .replace(".txt.gz", OUT_ENDING)
        .replace(".fastq", OUT_ENDING);
    let outdir = std::path::Path::new(&cfg.outdir);
    let outfile = outdir.join(oname);
    println!("{}", outfile.to_str().unwrap());

    let (reader, _) = niffler::from_path(infile).unwrap_formatful("Invalid input path!");
    let fq_reader = fastq::Reader::new(reader);
    for record in fq_reader.records() {
        let record = record.expect("Invalid fastq");
        println!("{}", record.id());
    }
}

struct ProgConfig {
    infile: String,
    outdir: String,
    rr_required: bool,
    qc_required: bool,
    left_flank: String,
    right_flank: String,
    insert_length: i32,
    expected_start: i32,
    position_tolerance: i32,
    min_peak_qual: u8,
    min_mean_qual: u8,
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
            min_peak_qual: json_config.qual_peak.unwrap_or(0), // TODO: Maybe change the config type to accept options
            min_mean_qual: json_config.qual_mean.unwrap_or(0),
        })
    }
}
