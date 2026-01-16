//! Enhanced error handling and diagnostics for CFD operations
//!
//! This module provides comprehensive error types, diagnostic information,
//! and recovery strategies for CFD simulations.

use std::collections::HashMap;
use std::fmt;

/// Enhanced error type with detailed diagnostic information
#[derive(Debug, Clone)]
pub struct CfdError {
    /// Error category
    pub category: ErrorCategory,
    /// Specific error code
    pub code: String,
    /// Human-readable message
    pub message: String,
    /// Technical details
    pub details: Option<String>,
    /// Suggested recovery actions
    pub recovery_actions: Vec<RecoveryAction>,
    /// Error severity
    pub severity: ErrorSeverity,
    /// Source location
    pub source_location: Option<SourceLocation>,
    /// Context information
    pub context: HashMap<String, String>,
}

/// Error categories for better error handling
#[derive(Debug, Clone, PartialEq)]
pub enum ErrorCategory {
    /// Input validation errors
    InputValidation,
    /// Numerical computation errors
    Numerical,
    /// Memory allocation errors
    Memory,
    /// Convergence errors
    Convergence,
    /// Mesh/geometry errors
    Geometry,
    /// File I/O errors
    Io,
    /// Parallel computation errors
    Parallel,
    /// GPU computation errors
    Gpu,
    /// Configuration errors
    Configuration,
    /// Physical constraint violations
    Physics,
    /// Performance degradation
    Performance,
    /// Unknown error type
    Unknown,
}

/// Error severity levels
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum ErrorSeverity {
    /// Informational message
    Info,
    /// Warning - simulation can continue
    Warning,
    /// Error - simulation should stop
    Error,
    /// Critical - immediate termination required
    Critical,
}

/// Recovery actions for error handling
#[derive(Debug, Clone)]
pub enum RecoveryAction {
    /// Retry the operation
    Retry { max_attempts: u32 },
    /// Reduce time step
    ReduceTimeStep { factor: f64 },
    /// Relax convergence criteria
    RelaxTolerance { factor: f64 },
    /// Switch to different solver
    SwitchSolver { solver_name: String },
    /// Refine mesh
    RefineMesh { factor: f64 },
    /// Check input parameters
    CheckInput,
    /// Restart from checkpoint
    RestartFromCheckpoint,
    /// Manual intervention required
    ManualIntervention { description: String },
}

/// Source location for error tracking
#[derive(Debug, Clone)]
pub struct SourceLocation {
    /// File name
    pub file: String,
    /// Line number
    pub line: u32,
    /// Function name
    pub function: String,
    /// Module path
    pub module: String,
}

impl CfdError {
    /// Create new error with basic information
    pub fn new(category: ErrorCategory, code: &str, message: &str) -> Self {
        Self {
            category,
            code: code.to_string(),
            message: message.to_string(),
            details: None,
            recovery_actions: Vec::new(),
            severity: ErrorSeverity::Error,
            source_location: None,
            context: HashMap::new(),
        }
    }
    
    /// Add technical details
    pub fn with_details(mut self, details: &str) -> Self {
        self.details = Some(details.to_string());
        self
    }
    
    /// Add recovery action
    pub fn with_recovery(mut self, action: RecoveryAction) -> Self {
        self.recovery_actions.push(action);
        self
    }
    
    /// Set severity
    pub fn with_severity(mut self, severity: ErrorSeverity) -> Self {
        self.severity = severity;
        self
    }
    
    /// Add source location
    pub fn with_location(mut self, file: &str, line: u32, function: &str, module: &str) -> Self {
        self.source_location = Some(SourceLocation {
            file: file.to_string(),
            line,
            function: function.to_string(),
            module: module.to_string(),
        });
        self
    }
    
    /// Add context information
    pub fn with_context(mut self, key: &str, value: &str) -> Self {
        self.context.insert(key.to_string(), value.to_string());
        self
    }
    
    /// Check if error is recoverable
    pub fn is_recoverable(&self) -> bool {
        !self.recovery_actions.is_empty()
    }
    
    /// Get recommended recovery action
    pub fn recommended_action(&self) -> Option<&RecoveryAction> {
        self.recovery_actions.first()
    }
}

impl fmt::Display for CfdError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{:?}] {}: {}", self.category, self.code, self.message)?;
        
        if let Some(details) = &self.details {
            write!(f, "\nDetails: {}", details)?;
        }
        
        if let Some(location) = &self.source_location {
            write!(f, "\nLocation: {}:{} in {}::{}", 
                   location.file, location.line, location.module, location.function)?;
        }
        
        if !self.recovery_actions.is_empty() {
            write!(f, "\nRecovery actions:")?;
            for (i, action) in self.recovery_actions.iter().enumerate() {
                write!(f, "\n  {}. {:?}", i + 1, action)?;
            }
        }
        
        Ok(())
    }
}

impl std::error::Error for CfdError {}

/// Convenience constructors for common error types
impl CfdError {
    /// Input validation error
    pub fn input_validation(message: &str) -> Self {
        Self::new(ErrorCategory::InputValidation, "INPUT_001", message)
            .with_recovery(RecoveryAction::CheckInput)
    }
    
    /// Numerical convergence error
    pub fn convergence_failure(iterations: u32, residual: f64) -> Self {
        Self::new(ErrorCategory::Convergence, "CONV_001", 
                 &format!("Failed to converge after {} iterations, residual: {:.2e}", 
                         iterations, residual))
            .with_details(format!("Final residual: {:.2e}, iterations: {}", residual, iterations).as_str())
            .with_recovery(RecoveryAction::RelaxTolerance { factor: 1.5 })
            .with_recovery(RecoveryAction::ReduceTimeStep { factor: 0.5 })
            .with_context("iterations", &iterations.to_string())
            .with_context("residual", &residual.to_string())
    }
    
    /// Memory allocation error
    pub fn memory_allocation(requested_bytes: usize) -> Self {
        Self::new(ErrorCategory::Memory, "MEM_001", 
                 &format!("Failed to allocate {} bytes", requested_bytes))
            .with_severity(ErrorSeverity::Critical)
            .with_context("requested_bytes", &requested_bytes.to_string())
    }
    
    /// Mesh quality error
    pub fn mesh_quality(issue: &str, cell_id: usize) -> Self {
        Self::new(ErrorCategory::Geometry, "MESH_001", 
                 &format!("Mesh quality issue in cell {}: {}", cell_id, issue))
            .with_recovery(RecoveryAction::RefineMesh { factor: 2.0 })
            .with_context("cell_id", &cell_id.to_string())
    }
    
    /// GPU computation error
    pub fn gpu_computation(operation: &str, device: &str) -> Self {
        Self::new(ErrorCategory::Gpu, "GPU_001", 
                 &format!("GPU computation failed for {} on device {}", operation, device))
            .with_recovery(RecoveryAction::SwitchSolver { 
                solver_name: "CPU".to_string() 
            })
            .with_context("operation", operation)
            .with_context("device", device)
    }
    
    /// Performance degradation warning
    pub fn performance_degradation(operation: &str, expected_ms: f64, actual_ms: f64) -> Self {
        let slowdown = actual_ms / expected_ms;
        Self::new(ErrorCategory::Performance, "PERF_001", 
                 &format!("Performance degradation in {}: {:.1}x slower than expected", 
                         operation, slowdown))
            .with_severity(ErrorSeverity::Warning)
            .with_details(format!("Expected: {:.2}ms, Actual: {:.2}ms, Slowdown: {:.1}x", 
                               expected_ms, actual_ms, slowdown).as_str())
            .with_context("operation", operation)
            .with_context("slowdown", &slowdown.to_string())
    }
    
    /// Physical constraint violation
    pub fn physics_violation(constraint: &str, value: f64, expected_range: &str) -> Self {
        Self::new(ErrorCategory::Physics, "PHYS_001", 
                 &format!("Physical constraint violation: {} = {} (expected: {})", 
                         constraint, value, expected_range))
            .with_severity(ErrorSeverity::Error)
            .with_context("constraint", constraint)
            .with_context("value", &value.to_string())
            .with_context("expected_range", expected_range)
    }
}

/// Error handler with automatic recovery strategies
pub struct ErrorHandler {
    /// Error history
    error_history: Vec<CfdError>,
    /// Recovery attempt limits
    recovery_limits: HashMap<String, u32>,
    /// Current recovery attempts
    recovery_attempts: HashMap<String, u32>,
}

impl ErrorHandler {
    /// Create new error handler
    pub fn new() -> Self {
        Self {
            error_history: Vec::new(),
            recovery_limits: HashMap::new(),
            recovery_attempts: HashMap::new(),
        }
    }
    
    /// Handle an error with automatic recovery
    pub fn handle_error(&mut self, error: CfdError) -> Result<(), CfdError> {
        self.error_history.push(error.clone());
        
        // Check if we've exceeded recovery attempts
        let error_key = format!("{:?}:{}", error.category, error.code);
        let attempts = self.recovery_attempts.entry(error_key.clone()).or_insert(0);
        let limit = self.recovery_limits.get(&error_key).copied().unwrap_or(3);
        
        if *attempts >= limit {
            return Err(error.with_severity(ErrorSeverity::Critical)
                          .with_details("Maximum recovery attempts exceeded"));
        }
        
        // Attempt recovery
        if let Some(action) = error.recommended_action() {
            *attempts += 1;
            println!("Attempting recovery: {:?}", action);
            // TODO: Implement actual recovery logic
            return Ok(());
        }
        
        Err(error)
    }
    
    /// Get error statistics
    pub fn get_error_statistics(&self) -> ErrorStatistics {
        let mut category_counts = HashMap::new();
        let mut severity_counts = HashMap::new();
        
        for error in &self.error_history {
            *category_counts.entry(error.category.clone()).or_insert(0) += 1;
            *severity_counts.entry(error.severity.clone()).or_insert(0) += 1;
        }
        
        ErrorStatistics {
            total_errors: self.error_history.len(),
            category_counts,
            severity_counts,
            most_common_category: self.find_most_common(&category_counts),
        }
    }
    
    /// Find most common error category
    fn find_most_common(&self, counts: &HashMap<ErrorCategory, u32>) -> Option<ErrorCategory> {
        counts.iter()
            .max_by_key(|(_, &count)| count)
            .map(|(category, _)| category.clone())
    }
    
    /// Clear error history
    pub fn clear_history(&mut self) {
        self.error_history.clear();
        self.recovery_attempts.clear();
    }
}

impl Default for ErrorHandler {
    fn default() -> Self {
        Self::new()
    }
}

/// Error statistics summary
#[derive(Debug, Clone)]
pub struct ErrorStatistics {
    /// Total number of errors
    pub total_errors: usize,
    /// Counts by category
    pub category_counts: HashMap<ErrorCategory, u32>,
    /// Counts by severity
    pub severity_counts: HashMap<ErrorSeverity, u32>,
    /// Most common error category
    pub most_common_category: Option<ErrorCategory>,
}

/// Macro for creating errors with source location
#[macro_export]
macro_rules! cfd_error {
    ($category:expr, $code:expr, $message:expr) => {
        CfdError::new($category, $code, $message)
            .with_location(file!(), line!(), function_name!(), module_path!())
    };
    ($category:expr, $code:expr, $message:expr, $($key:ident = $value:expr),*) => {
        cfd_error!($category, $code, $message)
            $(.with_context(stringify!($key), &$value.to_string()))*
    };
}

/// Get function name (requires nightly Rust)
#[cfg(not(nightly))]
fn function_name() -> &'static str {
    "unknown"
}

#[cfg(nightly)]
fn function_name() -> &'static str {
    fn __inner() -> &'static str {
        std::any::type_name::<fn()>()
            .split("::")
            .nth_back(1)
            .unwrap_or("unknown")
    }
    __inner()
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_error_creation() {
        let error = CfdError::convergence_failure(100, 1e-6);
        assert_eq!(error.category, ErrorCategory::Convergence);
        assert_eq!(error.code, "CONV_001");
        assert!(error.is_recoverable());
    }
    
    #[test]
    fn test_error_handler() {
        let mut handler = ErrorHandler::new();
        let error = CfdError::input_validation("Invalid parameter");
        
        // First attempt should succeed
        assert!(handler.handle_error(error.clone()).is_ok());
        
        // After limit should fail
        for _ in 0..3 {
            let _ = handler.handle_error(error.clone());
        }
        assert!(handler.handle_error(error).is_err());
    }
    
    #[test]
    fn test_error_statistics() {
        let mut handler = ErrorHandler::new();
        
        handler.handle_error(CfdError::input_validation("Test 1")).unwrap();
        handler.handle_error(CfdError::input_validation("Test 2")).unwrap();
        handler.handle_error(CfdError::convergence_failure(50, 1e-5)).unwrap();
        
        let stats = handler.get_error_statistics();
        assert_eq!(stats.total_errors, 3);
        assert_eq!(stats.category_counts.get(&ErrorCategory::InputValidation), Some(&2));
        assert_eq!(stats.category_counts.get(&ErrorCategory::Convergence), Some(&1));
    }
}
