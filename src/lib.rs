//! Tool to deal with targeted amplicon sequencing results
pub mod config;
pub mod handling;
pub mod match_type;
pub mod matching;
pub mod output;
pub mod stat;
#[macro_use]
extern crate clap;

#[allow(unused_imports)]
use bio::pattern_matching::{bndm, shift_and};
pub type ExactPattern = shift_and::ShiftAnd;
