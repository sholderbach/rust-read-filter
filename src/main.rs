use std::io;
use bio::io::fastq;
#[macro_use]
extern crate clap;

fn main() {
    // Specify CLI
    let matches = clap_app!(readfilter =>
        (version: "0.1")
        (author: "Stefan Holderbach")
        (about: "Read filter for amplicon sequencing with a defined region")
        (@arg CONFIG: -c --config +takes_value +required "Sets a custom config file")
        (@arg INPUT: +required "Sets the input file to use")
        (@arg OUPUT: +required "Sets the output path")
        (@arg read_report: -r --("read-report") "Also output a table with QC information for each read")
        (@arg qc_report: -q --("qc-report") "Also output a table with overall QC information")
        (@arg debug: -d ... "Sets the level of debugging information")
    ).get_matches();

    // Unpack arguments
    println!("{:?}", matches);
    let infile = matches.value_of("INPUT");
    let outdir = matches.value_of("OUTPUT");
    let config_file = matches.value_of("CONFIG");
    let rr_required = matches.is_present("read_report");
    let qc_required = matches.is_present("qc_report");
    // Same as before...
    
}