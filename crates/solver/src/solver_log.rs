//! Solver logging infrastructure.
//!
//! Captures solver log messages and writes them to a log file alongside
//! solver results.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use crate::error::SolverError;

/// Thread-safe solver log that accumulates messages and can flush to disk.
pub struct SolverLog {
    entries: Mutex<Vec<LogEntry>>,
    log_path: PathBuf,
}

/// A single log entry with timestamp and message.
#[derive(Debug, Clone)]
pub struct LogEntry {
    pub timestamp_secs: u64,
    pub level: LogLevel,
    pub message: String,
}

#[derive(Debug, Clone, Copy)]
pub enum LogLevel {
    Info,
    Warning,
    Error,
}

impl SolverLog {
    /// Create a new solver log that will be written to `work_dir/solver.log`.
    pub fn new(work_dir: &Path) -> Self {
        SolverLog {
            entries: Mutex::new(Vec::new()),
            log_path: work_dir.join("solver.log"),
        }
    }

    /// Add an info message.
    pub fn info(&self, message: impl Into<String>) {
        self.add(LogLevel::Info, message.into());
    }

    /// Add a warning message.
    pub fn warn(&self, message: impl Into<String>) {
        self.add(LogLevel::Warning, message.into());
    }

    /// Add an error message.
    pub fn error(&self, message: impl Into<String>) {
        self.add(LogLevel::Error, message.into());
    }

    fn add(&self, level: LogLevel, message: String) {
        let timestamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default()
            .as_secs();

        if let Ok(mut entries) = self.entries.lock() {
            entries.push(LogEntry {
                timestamp_secs: timestamp,
                level,
                message,
            });
        }
    }

    /// Flush all log entries to the log file.
    pub fn flush(&self) -> Result<(), SolverError> {
        let entries = self
            .entries
            .lock()
            .map_err(|_| SolverError::SolverExecution("Failed to lock log".into()))?;

        if entries.is_empty() {
            return Ok(());
        }

        let mut file = std::fs::File::create(&self.log_path)
            .map_err(|e| SolverError::io(&self.log_path, e))?;

        for entry in entries.iter() {
            let level_str = match entry.level {
                LogLevel::Info => "INFO",
                LogLevel::Warning => "WARN",
                LogLevel::Error => "ERROR",
            };
            writeln!(
                file,
                "[{}] [{}] {}",
                entry.timestamp_secs, level_str, entry.message
            )
            .map_err(|e| SolverError::io(&self.log_path, e))?;
        }

        Ok(())
    }

    /// Get all entries as a formatted string.
    pub fn to_string(&self) -> String {
        let entries = match self.entries.lock() {
            Ok(e) => e,
            Err(_) => return String::new(),
        };

        entries
            .iter()
            .map(|e| {
                let level = match e.level {
                    LogLevel::Info => "INFO",
                    LogLevel::Warning => "WARN",
                    LogLevel::Error => "ERROR",
                };
                format!("[{}] {}", level, e.message)
            })
            .collect::<Vec<_>>()
            .join("\n")
    }

    /// Get the path to the log file.
    pub fn log_path(&self) -> &Path {
        &self.log_path
    }
}
