//! Utilities for loading user config from the command line and json config files
use crate::handling::{GracefulOption, GracefulResult};
use serde::Deserialize;
use std::{error::Error, fs::File};

/// JSON config for read_filter
/// ## Example
/// ``` json
/// {
///     "left_flank": "AGAGAGGC",
///     "right_flank": "GCCCAGGC",
///     "content_length": 21,
///     "expect_begin": 36,
///     "tolerance": 100,
///     "qual_peak": 20,
///     "qual_mean": 30
/// }
/// ```
#[derive(Deserialize)]
pub struct FilterConf {
    pub left_flank: String,
    pub right_flank: String,
    pub content_length: u32,
    pub expect_begin: u32,
    pub tolerance: u32,
    pub qual_peak: Option<u8>,
    pub qual_mean: Option<u8>,
}

pub fn load_json_config<P: AsRef<std::path::Path>>(
    json_path: P,
) -> Result<FilterConf, Box<dyn std::error::Error>> {
    let reader = File::open(json_path)?;
    let res: FilterConf = serde_json::from_reader(reader)?;
    Ok(res)
}
#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_json_type() {
        let js_str: &'static str = r#"{
      "left_flank": "AGGGCCAG",
      "right_flank": "GCCCAGGC",
      "content_length": 27,
      "expect_begin": 36,
      "tolerance": 8,
      "qual_peak":20,
      "qual_mean":30
  }"#;

        let result: FilterConf = serde_json::from_str(js_str).unwrap();
        assert_eq!(result.qual_mean.unwrap(), 30u8);
    }
}

/// Summarized config used by different parts of the program
pub struct ProgConfig {
    pub infile: String,
    pub outdir: String,
    pub rr_required: bool,
    pub qc_required: bool,
    pub left_flank: String,
    pub right_flank: String,
    pub insert_length: u32,
    pub expected_start: u32,
    pub position_tolerance: u32,
    pub min_peak_qual: Option<u8>,
    pub min_mean_qual: Option<u8>,
}

impl ProgConfig {
    pub fn from_cli() -> Result<ProgConfig, Box<dyn Error>> {
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
            load_json_config(config_file).unwrap_formatful("While parsing configuration");

        // TODO: Add checks to block useless or illegal inputs/configs
        // Category illegal:
        // pattern larger than 64 (shift_and or bndm wouldn't be able to construct the match vector)
        // Category useless:
        // If read length would be known, exptected start + total pattern length beyond read_length

        // TODO: Make positional limits optional (assumption only the valid entity will match the full pattern and quality based ranking is unnecessary)

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
