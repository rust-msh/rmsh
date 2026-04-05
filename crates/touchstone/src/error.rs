use thiserror::Error;

#[derive(Debug, Clone, Error)]
pub enum TouchstoneError {
    #[error("invalid option line: {0}")]
    InvalidOptionLine(String),

    #[error("line {line}: expected {expected} values, got {got}")]
    InvalidDataLine {
        line: usize,
        expected: usize,
        got: usize,
    },

    #[error("line {line}: invalid number \"{value}\"")]
    InvalidNumber { line: usize, value: String },

    #[error("line {line}: {message}")]
    ParseError { line: usize, message: String },

    #[error("unsupported version: {0}")]
    UnsupportedVersion(String),

    #[error("no option line found (expected line starting with '#')")]
    NoOptionLine,

    #[error("no data found")]
    NoData,

    #[error("inconsistent port count: header says {expected}, data has {got}")]
    InconsistentPortCount { expected: usize, got: usize },
}
