pub mod handling {
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
}
pub mod jsonconf {
    use serde::Deserialize;
    use std::fs::File;

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

    pub fn load_json_config(json_path: &str) -> Result<FilterConf, Box<dyn std::error::Error>> {
        let reader = File::open(json_path)?;
        let res: FilterConf = serde_json::from_reader(reader)?;
        Ok(res)
    }
}
#[cfg(test)]
mod test {
    use super::jsonconf;

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

        let result: jsonconf::FilterConf = serde_json::from_str(js_str).unwrap();
        assert_eq!(result.qual_mean.unwrap(), 30u8);
    }
}
