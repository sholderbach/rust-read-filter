use bio::io::fastq;
// TODO: Decide on json crate for config file

#[macro_use]
extern crate clap;

#[allow(unused)]
fn main() {
    // May be consts
    let out_ending = ".processed.tsv";
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
    let infile = matches.value_of("INPUT").unwrap();
    let outdir = matches.value_of("OUTPUT").unwrap();
    let config_file = matches.value_of("CONFIG").unwrap();
    let rr_required = matches.is_present("read_report");
    let qc_required = matches.is_present("qc_report");
    // Same as before...
    let infile = std::path::Path::new(infile);
    let bname = infile.file_name().expect("Input needs to be a file");
    let oname = String::from(bname.to_str().unwrap()).replace(".fastq.gz", out_ending)
                                                     .replace(".txt.gz", out_ending)
                                                     .replace(".fastq", out_ending);
    let outdir = std::path::Path::new(outdir);
    let outfile = outdir.join(oname);
    println!("{}", outfile.to_str().unwrap());

    let (reader, _) = niffler::from_path(infile).expect("Invalid input path!");
    let fq_reader = fastq::Reader::new(reader);
    for record in fq_reader.records() {
        let record = record.expect("Invalid fastq");
        println!("{}", record.id());
    }
}