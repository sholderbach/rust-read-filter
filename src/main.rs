use bio::io::fastq;
use counter::Counter;
use read_filter::config::ProgConfig;
use read_filter::handling::{GracefulOption, GracefulResult};
use read_filter::matching::ReadFilter;
use read_filter::output::{write_config_header, write_read_report_header, write_stats_header};
use read_filter::stat::{QualStats, RunningStats};
#[allow(unused_imports)]
use std::todo;
use std::{
    fs::File,
    io::{BufWriter, Write},
    iter::Iterator,
    path::{Path, PathBuf},
};

const OUT_ENDING: &str = ".processed.tsv";
const RR_ENDING: &str = ".readreport.tsv";
const QC_ENDING: &str = ".quality.tsv";
fn main() {
    let cfg = ProgConfig::from_cli().unwrap_graceful();
    let infile = Path::new(&cfg.infile);
    let outdir = Path::new(&cfg.outdir);
    std::fs::create_dir_all(outdir).unwrap_messageful(&format!(
        "Could not create output directory at: {:?}",
        outdir.to_str().unwrap()
    ));
    let bname = infile
        .file_name()
        .unwrap_graceful("Input needs to be a file"); // TODO file name is just basename only errs when given ..
    let get_outpath = |ending: &str| -> PathBuf {
        let oname = String::from(bname.to_str().unwrap())
            .replace(".fastq.gz", ending)
            .replace(".txt.gz", ending)
            .replace(".fastq", ending);
        outdir.join(oname)
    };
    fn create_file(path: PathBuf) -> BufWriter<File> {
        let f = std::fs::File::create(&path).unwrap_messageful(&format!(
            "Could not create output file at: {:?}",
            path.to_str().unwrap()
        ));
        BufWriter::new(f)
    }
    let outfile = get_outpath(OUT_ENDING);
    let rr_file = get_outpath(RR_ENDING);
    let qc_file = get_outpath(QC_ENDING);
    let mut ofile = create_file(outfile);

    // FASTQ parsing
    let (reader, _compression) = niffler::from_path(infile).unwrap_formatful("Invalid input path!");
    let fq_reader = fastq::Reader::new(reader);

    let mut stats = RunningStats::default();
    let mut qual_stats = QualStats::new();
    let rf = ReadFilter::new(fq_reader.records(), &cfg, &mut stats);

    let counter = match (cfg.rr_required, cfg.qc_required) {
        (false, false) => rf.map(|a| a.seq).collect::<Counter<Vec<u8>>>(),
        (true, false) => {
            let mut rr_file = create_file(rr_file);
            write_read_report_header(&mut rr_file, cfg.insert_length as usize)
                .unwrap_messageful("Error while writing output");
            rf.map(|a| {
                a.write_read_report_line(&mut rr_file).unwrap();
                a.seq
            })
            .collect::<Counter<Vec<u8>>>()
        }
        (false, true) => rf
            .map(|a| {
                qual_stats.append(&a);
                a.seq
            })
            .collect::<Counter<Vec<u8>>>(),
        (true, true) => {
            let mut rr_file = create_file(rr_file);
            write_read_report_header(&mut rr_file, cfg.insert_length as usize)
                .unwrap_messageful("Error while writing output");
            rf.map(|a| {
                qual_stats.append(&a);
                a.write_read_report_line(&mut rr_file).unwrap();
                a.seq
            })
            .collect::<Counter<Vec<u8>>>()
        }
    };

    write_config_header(&mut ofile, &cfg).unwrap_messageful("Error while writing output");
    write_stats_header(&mut ofile, &stats).unwrap();
    writeln!(ofile, "seq\treads").unwrap();
    for (seq, count) in counter.iter() {
        writeln!(ofile, "{}\t{}", std::str::from_utf8(seq).unwrap(), count).unwrap();
    }

    if cfg.qc_required {
        let mut qc_file = create_file(qc_file);
        qual_stats
            .write_to_buf(&mut qc_file, cfg.insert_length as usize)
            .unwrap_messageful("Error while writing output");
    }
}
