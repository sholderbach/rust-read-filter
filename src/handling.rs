//! Contains Traits to unwrap without panic but a deliberate process exit.
//! Behavior regarding teardown destructors might need to be considered as it calls std::process:exit
// TODO Reconsider this due to potential teardown issues.
use std::process;
/// Lets the program exit with exit code one based on Result<T,E> with desired messages to stderr
pub trait GracefulResult<T> {
    fn unwrap_graceful(self) -> T;
    fn unwrap_messageful(self, msg: &str) -> T;
    fn unwrap_formatful(self, fmt: &str) -> T;
}

impl<T, E: std::fmt::Display> GracefulResult<T> for Result<T, E> {
    /// If Err(e) print e to stderr and exit with 1.
    fn unwrap_graceful(self) -> T {
        match self {
            Err(e) => {
                eprintln!("{}", e);
                process::exit(1);
            }
            Ok(x) => x,
        }
    }
    /// If Err(_) print msg to stderr and exit with 1.
    fn unwrap_messageful(self, msg: &str) -> T {
        match self {
            Err(_) => {
                eprintln!("{}", msg);
                process::exit(1);
            }
            Ok(x) => x,
        }
    }
    /// If Err(e) print "{msg}: {e}" to stderr and exit with 1.
    fn unwrap_formatful(self, fmt: &str) -> T {
        match self {
            Err(e) => {
                eprintln!("{}: {}", fmt, e);
                process::exit(1);
            }
            Ok(x) => x,
        }
    }
}

/// Shortcut to fail with exit code 1 and a custom message for empty options
pub trait GracefulOption<T> {
    fn unwrap_graceful(self, msg: &str) -> T;
}

impl<T> GracefulOption<T> for Option<T> {
    /// If None print msg to stderr and exit with 1.
    fn unwrap_graceful(self, msg: &str) -> T {
        match self {
            None => {
                eprintln!("{}", msg);
                process::exit(1);
            }
            Some(x) => x,
        }
    }
}
