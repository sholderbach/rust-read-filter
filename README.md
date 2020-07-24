# WIP: A rust amplicon sequencing "read_filter" tool

I wanted to learn rust with a useful bioinformatics project, and thus created this crate to reimplement an existing python tool of mine.
**So far the functionality is not implemented and this repo serves more to document my learning process.**
May squash the history of this repo later.

## Goal of the tool

The existing tool was able to filter sequence fragments from amplicon sequencing experiments based on exact match sequence context around an insert.
Additional quality filtering could be performed. The gathered fragments were counted and exported in a custom table for use in the analysis of directed evolution experiments.
Additionally the original tool was able to report on the quality distribution of the matching reads in a binned data format for a quality dashboard. For debug purposes my original tool was also able to export matching reads with associated qualities in a custom format. (Standard fastq might be a more civil choice for the future)

## Existing CLI

The CLI currently matches the existing one as the tool is part of snakemake pipelines. Might refactor to a more UNIXy layout at some point.
For easier integration into the pipeline, the filtering details are currently provided by a config json file.

```
USAGE:
    read_filter [FLAGS] <INPUT> <OUTPUT> --config <CONFIG>

FLAGS:
    -d                   Sets the level of debugging information
    -h, --help           Prints help information
    -q, --qc-report      Also output a table with overall QC information
    -r, --read-report    Also output a table with QC information for each read
    -V, --version        Prints version information

OPTIONS:
    -c, --config <CONFIG>    Sets a custom config file

ARGS:
    <INPUT>     Sets the input file to use
    <OUTPUT>    Sets the output path
```
