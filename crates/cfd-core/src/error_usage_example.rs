//! Example demonstrating proper usage of the improved error system

use crate::error_fixed::{*, Result};

/// Example: Loading and validating a solver configuration
pub fn load_solver_config(path: &str) -> Result<SolverConfig> {
    // Read the file - specific I/O error with context
    let contents = std::fs::read_to_string(path)
        .map_err(|source| IOError::ReadError {
            path: path.to_string(),
            source,
        })?;
    
    // Parse JSON - structured parsing error
    let config: serde_json::Value = serde_json::from_str(&contents)
        .map_err(|e| IOError::ParseError {
            path: path.to_string(),
            line: e.line(),
            source: Box::new(e),
        })?;
    
    // Validate required fields - structured configuration errors
    let tolerance = config["tolerance"]
        .as_f64()
        .ok_or_else(|| ConfigurationError::MissingField {
            field: "tolerance".into(),
        })?;
    
    if tolerance <= 0.0 {
        return Err(ConfigurationError::InvalidValue {
            field: "tolerance".into(),
            reason: "must be positive".into(),
        }.into());
    }
    
    let max_iterations = config["max_iterations"]
        .as_u64()
        .ok_or_else(|| ConfigurationError::MissingField {
            field: "max_iterations".into(),
        })? as usize;
    
    if max_iterations == 0 {
        return Err(ConfigurationError::InvalidValue {
            field: "max_iterations".into(),
            reason: "must be at least 1".into(),
        }.into());
    }
    
    Ok(SolverConfig {
        tolerance,
        max_iterations,
    })
}

/// Example: Matrix operations with proper error handling
pub fn solve_linear_system(
    matrix: &Matrix,
    rhs: &Vector,
) -> Result<Vector> {
    // Check dimensions - structured input error
    if matrix.cols() != rhs.len() {
        return Err(InputError::DimensionMismatch {
            expected: matrix.cols(),
            actual: rhs.len(),
        }.into());
    }
    
    // Check for singular matrix - specific numerical error
    let cond = matrix.condition_number();
    if cond > 1e15 {
        return Err(NumericalError::SingularMatrix { cond }.into());
    }
    
    // Attempt to solve
    let mut solver = LinearSolver::new();
    let solution = solver.solve(matrix, rhs)
        .map_err(|_| LinearSolverError::ConvergenceFailed {
            iterations: solver.iteration_count(),
        })?;
    
    Ok(solution)
}

/// Example: Iterative solver with recoverable error handling
pub fn iterative_solve_with_recovery(
    problem: &Problem,
    initial_guess: &Solution,
) -> Result<Solution> {
    let mut solution = initial_guess.clone();
    let mut max_iters = 100;
    let mut attempts = 0;
    const MAX_ATTEMPTS: usize = 3;
    
    loop {
        match iterative_solve_internal(problem, &solution, max_iters) {
            Ok(converged) => return Ok(converged),
            
            Err(Error::Convergence(e)) if e.is_recoverable() => {
                // Handle recoverable convergence errors
                attempts += 1;
                if attempts >= MAX_ATTEMPTS {
                    return Err(e.into());
                }
                
                match e {
                    ConvergenceError::MaxIterationsExceeded { .. } => {
                        // Double the iteration limit
                        max_iters *= 2;
                        println!("Increasing max iterations to {}", max_iters);
                    }
                    ConvergenceError::Stagnation { .. } => {
                        // Perturb the solution slightly
                        solution.perturb(1e-6);
                        println!("Perturbing solution to escape stagnation");
                    }
                    _ => return Err(e.into()),
                }
            }
            
            Err(e) => return Err(e),
        }
    }
}

/// Example: Using the error system for validation
pub fn validate_mesh(mesh: &Mesh) -> Result<()> {
    // Check for degenerate elements
    for (i, element) in mesh.elements().enumerate() {
        let quality = element.quality_metric();
        if quality < 1e-10 {
            return Err(MeshError::DegenerateElement {
                index: i,
                reason: format!("quality metric {:.3e} below threshold", quality),
            }.into());
        }
    }
    
    // Check overall mesh quality
    let avg_quality = mesh.average_quality();
    const QUALITY_THRESHOLD: f64 = 0.3;
    if avg_quality < QUALITY_THRESHOLD {
        return Err(MeshError::PoorQuality {
            metric: "average_quality".into(),
            value: avg_quality,
            threshold: QUALITY_THRESHOLD,
        }.into());
    }
    
    Ok(())
}

/// Example: Programmatic error handling by consumers
pub fn consumer_example() {
    let result = load_solver_config("config.json");
    
    match result {
        Ok(config) => {
            println!("Config loaded successfully");
        }
        Err(Error::Configuration(ConfigurationError::MissingField { field })) => {
            // Specific handling for missing fields
            eprintln!("Please add '{}' to your configuration file", field);
        }
        Err(Error::Configuration(ConfigurationError::InvalidValue { field, reason })) => {
            // Specific handling for invalid values
            eprintln!("Invalid value for '{}': {}", field, reason);
        }
        Err(Error::IO(IOError::ReadError { path, .. })) => {
            // Specific handling for file not found
            eprintln!("Could not read configuration file: {}", path);
        }
        Err(e) => {
            // Generic handling for other errors
            eprintln!("Error: {}", e);
            
            // Check if we can provide recovery hints
            if let Some(hint) = e.recovery_hint() {
                eprintln!("Hint: {}", hint);
            }
        }
    }
}

// Placeholder types for the example
struct SolverConfig {
    tolerance: f64,
    max_iterations: usize,
}

struct Matrix;
impl Matrix {
    fn cols(&self) -> usize { 10 }
    fn condition_number(&self) -> f64 { 1.0 }
}

struct Vector;
impl Vector {
    fn len(&self) -> usize { 10 }
}

struct LinearSolver {
    iterations: usize,
}

impl LinearSolver {
    fn new() -> Self { Self { iterations: 0 } }
    fn solve(&mut self, _: &Matrix, _: &Vector) -> std::result::Result<Vector, ()> {
        self.iterations = 50;
        Ok(Vector)
    }
    fn iteration_count(&self) -> usize { self.iterations }
}

struct Problem;
struct Solution;
impl Solution {
    fn clone(&self) -> Self { Solution }
    fn perturb(&mut self, _: f64) {}
}

fn iterative_solve_internal(_: &Problem, _: &Solution, _: usize) -> Result<Solution> {
    Err(ConvergenceError::MaxIterationsExceeded {
        max: 100,
        residual: 1e-3,
    }.into())
}

struct Mesh;
impl Mesh {
    fn elements(&self) -> impl Iterator<Item = Element> {
        vec![Element].into_iter()
    }
    fn average_quality(&self) -> f64 { 0.5 }
}

struct Element;
impl Element {
    fn quality_metric(&self) -> f64 { 0.5 }
}