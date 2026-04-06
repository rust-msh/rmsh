use thiserror::Error;

#[derive(Debug, Error)]
pub enum AnsysExchangeError {
    #[error("IO error: {0}")]
    Io(String),
    #[error("parse error: {0}")]
    Parse(String),
    #[error("unsupported content: {0}")]
    Unsupported(String),
    #[error("JSON error: {0}")]
    Json(String),
}
