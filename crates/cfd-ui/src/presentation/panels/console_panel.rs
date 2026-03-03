//! Console panel — log output and solver messages.

/// A log message entry.
#[derive(Clone, Debug)]
pub struct LogEntry {
    /// Log level.
    pub level: LogLevel,
    /// Message text.
    pub message: String,
}

/// Log severity level.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LogLevel {
    Info,
    Warning,
    Error,
}

/// Console state holding log messages.
#[derive(Clone, Debug, Default)]
pub struct ConsoleState {
    entries: Vec<LogEntry>,
    max_entries: usize,
}

impl ConsoleState {
    /// Create a new console with a maximum entry limit.
    #[must_use]
    pub fn new(max_entries: usize) -> Self {
        Self {
            entries: Vec::new(),
            max_entries,
        }
    }

    /// Add a log message.
    pub fn log(&mut self, level: LogLevel, message: String) {
        self.entries.push(LogEntry { level, message });
        if self.entries.len() > self.max_entries {
            self.entries.remove(0);
        }
    }

    /// Add an info message.
    pub fn info(&mut self, message: String) {
        self.log(LogLevel::Info, message);
    }

    /// Add a warning message.
    pub fn warn(&mut self, message: String) {
        self.log(LogLevel::Warning, message);
    }

    /// Add an error message.
    pub fn error(&mut self, message: String) {
        self.log(LogLevel::Error, message);
    }

    /// All log entries.
    #[must_use]
    pub fn entries(&self) -> &[LogEntry] {
        &self.entries
    }

    /// Clear all entries.
    pub fn clear(&mut self) {
        self.entries.clear();
    }
}
