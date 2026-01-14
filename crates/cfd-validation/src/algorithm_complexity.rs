//! Comprehensive algorithm complexity analysis for CFD operations
//!
//! This module provides detailed Big-O complexity analysis and performance
//! characterization for all major CFD algorithms implemented in the codebase.
//!
//! ## Algorithm Complexity Registry
//!
//! This registry documents the theoretical and practical performance characteristics
//! of key CFD algorithms, enabling informed algorithm selection and optimization.

use std::collections::HashMap;
use std::fmt::Write as _;

/// Algorithm complexity information
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct AlgorithmComplexityInfo {
    /// Algorithm name
    pub name: String,
    /// Time complexity (Big-O notation)
    pub time_complexity: String,
    /// Space complexity (Big-O notation)
    pub space_complexity: String,
    /// Memory access pattern description
    pub memory_pattern: String,
    /// Cache efficiency rating (0.0 to 1.0)
    pub cache_efficiency: f64,
    /// Parallel scalability factor (0.0 to 1.0)
    pub scalability: f64,
    /// Literature references for complexity analysis
    pub references: Vec<String>,
    /// Performance notes and optimization hints
    pub notes: Vec<String>,
}

impl AlgorithmComplexityInfo {
    /// Create new algorithm complexity info
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            time_complexity: "O(?)".to_string(),
            space_complexity: "O(?)".to_string(),
            memory_pattern: "Unknown".to_string(),
            cache_efficiency: 0.5,
            scalability: 0.5,
            references: Vec::new(),
            notes: Vec::new(),
        }
    }

    /// Set time complexity
    pub fn time_complexity(mut self, complexity: &str) -> Self {
        self.time_complexity = complexity.to_string();
        self
    }

    /// Set space complexity
    pub fn space_complexity(mut self, complexity: &str) -> Self {
        self.space_complexity = complexity.to_string();
        self
    }

    /// Set memory access pattern
    pub fn memory_pattern(mut self, pattern: &str) -> Self {
        self.memory_pattern = pattern.to_string();
        self
    }

    /// Set cache efficiency (0.0 = poor, 1.0 = excellent)
    pub fn cache_efficiency(mut self, efficiency: f64) -> Self {
        self.cache_efficiency = efficiency.clamp(0.0, 1.0);
        self
    }

    /// Set parallel scalability (0.0 = no speedup, 1.0 = perfect scaling)
    pub fn scalability(mut self, scalability: f64) -> Self {
        self.scalability = scalability.clamp(0.0, 1.0);
        self
    }

    /// Add literature reference
    pub fn reference(mut self, reference: &str) -> Self {
        self.references.push(reference.to_string());
        self
    }

    /// Add performance note
    pub fn note(mut self, note: &str) -> Self {
        self.notes.push(note.to_string());
        self
    }
}

/// Comprehensive CFD algorithm complexity registry
pub struct AlgorithmComplexityRegistry {
    algorithms: HashMap<String, AlgorithmComplexityInfo>,
}

impl AlgorithmComplexityRegistry {
    /// Create the complete CFD algorithm complexity registry
    pub fn new() -> Self {
        let mut registry = HashMap::new();

        // Linear solvers
        registry.insert(
            "ConjugateGradient".to_string(),
            AlgorithmComplexityInfo::new("ConjugateGradient")
                .time_complexity("O(N^{3/2})")
                .space_complexity("O(N²)")
                .memory_pattern("Sparse matrix-vector products with irregular access")
                .cache_efficiency(0.7)
                .scalability(0.8)
                .reference("Saad (2003): Iterative Methods for Sparse Linear Systems")
                .reference("Barrett et al. (1994): Templates for the Solution of Linear Systems")
                .note("SPMV bottleneck: memory bandwidth critical")
                .note("Preconditioning reduces iterations from O(N) to O(√N)"),
        );

        registry.insert(
            "GMRES".to_string(),
            AlgorithmComplexityInfo::new("GMRES")
                .time_complexity("O(N^{3/2})")
                .space_complexity("O(N²)")
                .memory_pattern("Arnoldi orthogonalization with dense matrix operations")
                .cache_efficiency(0.6)
                .scalability(0.75)
                .reference("Saad & Schultz (1986): GMRES: A generalized minimal residual algorithm")
                .note("Restarted GMRES (GMRES(m)) reduces memory to O(m*N)")
                .note("Orthogonalization dominates compute time"),
        );

        registry.insert(
            "BiCGSTAB".to_string(),
            AlgorithmComplexityInfo::new("BiCGSTAB")
                .time_complexity("O(N^{3/2})")
                .space_complexity("O(N²)")
                .memory_pattern(
                    "Bi-orthogonalization with two matrix-vector products per iteration",
                )
                .cache_efficiency(0.65)
                .scalability(0.8)
                .reference(
                    "van der Vorst (1992): Bi-CGSTAB: A fast and smoothly converging variant",
                )
                .note("More stable than CGS for non-symmetric systems")
                .note("Two SPMV operations per iteration"),
        );

        registry.insert(
            "Multigrid".to_string(),
            AlgorithmComplexityInfo::new("Multigrid")
                .time_complexity("O(N)")
                .space_complexity("O(N)")
                .memory_pattern("Hierarchical grid operations with smoothing and restriction")
                .cache_efficiency(0.75)
                .scalability(0.7)
                .reference("Trottenberg et al. (2001): Multigrid methods")
                .reference("Briggs et al. (2000): A multigrid tutorial")
                .note("Optimal complexity: O(N) vs O(N^{3/2}) for Krylov methods")
                .note("Setup cost amortized over many solves"),
        );

        // Time integration schemes
        registry.insert(
            "RungeKutta4".to_string(),
            AlgorithmComplexityInfo::new("RungeKutta4")
                .time_complexity("O(N)")
                .space_complexity("O(N)")
                .memory_pattern("Sequential stage computations with vector operations")
                .cache_efficiency(0.85)
                .scalability(0.9)
                .reference("Hairer & Nørsett (1993): Solving Ordinary Differential Equations I")
                .note("4 RHS evaluations per step")
                .note("CFL limit ≈ 2.8 for linear advection"),
        );

        registry.insert(
            "RungeKutta3".to_string(),
            AlgorithmComplexityInfo::new("RungeKutta3")
                .time_complexity("O(N)")
                .space_complexity("O(N)")
                .memory_pattern("Sequential stage computations with vector operations")
                .cache_efficiency(0.85)
                .scalability(0.9)
                .reference("Kennedy & Carpenter (2003): Additive Runge-Kutta schemes")
                .note("3 RHS evaluations per step")
                .note("CFL limit ≈ 1.7, good for nonlinear problems"),
        );

        registry.insert(
            "LowStorageRK4".to_string(),
            AlgorithmComplexityInfo::new("LowStorageRK4")
                .time_complexity("O(N)")
                .space_complexity("O(N)")
                .memory_pattern("In-place stage updates with minimal memory overhead")
                .cache_efficiency(0.9)
                .scalability(0.95)
                .reference("Bijl & Carpenter (2009): Low-order Runge-Kutta methods for CFD")
                .note("Memory efficient: O(N) vs O(5N) for classical RK4")
                .note("Same CFL limit as classical RK4"),
        );

        // Turbulence models
        registry.insert(
            "KEpsilon".to_string(),
            AlgorithmComplexityInfo::new("KEpsilon")
                .time_complexity("O(N)")
                .space_complexity("O(N)")
                .memory_pattern("Point-wise operations on turbulence variables")
                .cache_efficiency(0.9)
                .scalability(0.95)
                .reference(
                    "Launder & Spalding (1974): The numerical computation of turbulent flows",
                )
                .note("Two transport equations: k and ε")
                .note("Gradient computations dominate cost"),
        );

        registry.insert(
            "KOmegaSST".to_string(),
            AlgorithmComplexityInfo::new("KOmegaSST")
                .time_complexity("O(N)")
                .space_complexity("O(N)")
                .memory_pattern("Point-wise operations with cross-diffusion terms")
                .cache_efficiency(0.9)
                .scalability(0.95)
                .reference("Menter (1994): Two-equation eddy-viscosity turbulence models")
                .note("Superior near-wall behavior vs k-ε")
                .note("Additional cross-diffusion term increases complexity"),
        );

        registry.insert(
            "SpalartAllmaras".to_string(),
            AlgorithmComplexityInfo::new("SpalartAllmaras")
                .time_complexity("O(N)")
                .space_complexity("O(N)")
                .memory_pattern("Single transport equation with wall distance")
                .cache_efficiency(0.9)
                .scalability(0.95)
                .reference("Spalart & Allmaras (1994): A one-equation turbulence model")
                .note("One-equation model, simpler than two-equation")
                .note("Excellent for aerospace applications"),
        );

        registry.insert(
            "WALE".to_string(),
            AlgorithmComplexityInfo::new("WALE")
                .time_complexity("O(N)")
                .space_complexity("O(N)")
                .memory_pattern("Strain rate tensor computation with velocity gradients")
                .cache_efficiency(0.85)
                .scalability(0.9)
                .reference("Nicoud & Ducros (1999): Subgrid-scale stress modelling")
                .note("Superior near-wall behavior for LES")
                .note("More expensive than Smagorinsky due to tensor operations"),
        );

        registry.insert(
            "Smagorinsky".to_string(),
            AlgorithmComplexityInfo::new("Smagorinsky")
                .time_complexity("O(N)")
                .space_complexity("O(N)")
                .memory_pattern("Strain rate magnitude computation")
                .cache_efficiency(0.9)
                .scalability(0.95)
                .reference("Smagorinsky (1963): General circulation experiments")
                .note("Simple and efficient eddy viscosity model")
                .note("Isotropic assumption may be too restrictive"),
        );

        // Spatial discretization
        registry.insert(
            "FiniteDifference".to_string(),
            AlgorithmComplexityInfo::new("FiniteDifference")
                .time_complexity("O(N)")
                .space_complexity("O(N)")
                .memory_pattern("Stencil operations on structured grids")
                .cache_efficiency(0.95)
                .scalability(0.98)
                .reference("Hirsch (2007): Numerical computation of internal and external flows")
                .note("Excellent cache performance on structured grids")
                .note("Limited to simple geometries"),
        );

        registry.insert(
            "FiniteVolume".to_string(),
            AlgorithmComplexityInfo::new("FiniteVolume")
                .time_complexity("O(N)")
                .space_complexity("O(N)")
                .memory_pattern("Face-based flux computations")
                .cache_efficiency(0.8)
                .scalability(0.9)
                .reference("LeVeque (2002): Finite Volume Methods for Hyperbolic Problems")
                .note("Conservative by construction")
                .note("Flexible for complex geometries"),
        );

        registry.insert(
            "DiscontinuousGalerkin".to_string(),
            AlgorithmComplexityInfo::new("DiscontinuousGalerkin")
                .time_complexity("O(N * p^d)")
                .space_complexity("O(N * p^d)")
                .memory_pattern("Element-local operations with surface flux exchanges")
                .cache_efficiency(0.7)
                .scalability(0.8)
                .reference("Cockburn & Shu (2001): Runge-Kutta discontinuous Galerkin methods")
                .note("Arbitrary order accuracy with p refinement")
                .note("Higher memory per degree of freedom"),
        );

        // Validation and analysis
        registry.insert(
            "RichardsonExtrapolation".to_string(),
            AlgorithmComplexityInfo::new("RichardsonExtrapolation")
                .time_complexity("O(N * M)")
                .space_complexity("O(N)")
                .memory_pattern("Multiple grid evaluations with convergence analysis")
                .cache_efficiency(0.7)
                .scalability(0.8)
                .reference("Roache (1998): Verification and Validation in Computational Science")
                .note("M = number of grid levels for extrapolation")
                .note("Provides error estimates and convergence rates"),
        );

        registry.insert("ManufacturedSolutions".to_string(), AlgorithmComplexityInfo::new("ManufacturedSolutions")
            .time_complexity("O(N)")
            .space_complexity("O(1)")
            .memory_pattern("Analytical function evaluation at grid points")
            .cache_efficiency(0.95)
            .scalability(0.99)
            .reference("Salari & Knupp (2000): Code verification by the method of manufactured solutions")
            .note("Exact solution known for error analysis")
            .note("Function evaluation is typically cheap"));

        Self {
            algorithms: registry,
        }
    }

    /// Get complexity information for an algorithm
    pub fn get(&self, name: &str) -> Option<&AlgorithmComplexityInfo> {
        self.algorithms.get(name)
    }

    /// List all algorithms in the registry
    pub fn algorithms(&self) -> Vec<&AlgorithmComplexityInfo> {
        self.algorithms.values().collect()
    }

    /// Find algorithms by time complexity class
    pub fn by_time_complexity(&self, complexity: &str) -> Vec<&AlgorithmComplexityInfo> {
        self.algorithms
            .values()
            .filter(|algo| algo.time_complexity == complexity)
            .collect()
    }

    /// Find best algorithms for a given problem size
    pub fn recommend_for_size(
        &self,
        problem_size: usize,
    ) -> Vec<(&AlgorithmComplexityInfo, String)> {
        let mut recommendations = Vec::new();

        // For small problems (N < 10^4), higher-order methods are feasible
        if problem_size < 10_000 {
            if let Some(algo) = self.get("RungeKutta4") {
                recommendations.push((algo, "RK4 suitable for small problems".to_string()));
            }
            if let Some(algo) = self.get("DiscontinuousGalerkin") {
                recommendations.push((algo, "DG methods feasible for small N".to_string()));
            }
        }

        // For medium problems (10^4 < N < 10^7), efficient iterative methods
        if (10_000..10_000_000).contains(&problem_size) {
            if let Some(algo) = self.get("Multigrid") {
                recommendations.push((algo, "Multigrid optimal for medium problems".to_string()));
            }
            if let Some(algo) = self.get("ConjugateGradient") {
                recommendations.push((algo, "CG efficient with good preconditioning".to_string()));
            }
        }

        // For large problems (N > 10^7), scalability is critical
        if problem_size >= 10_000_000 {
            if let Some(algo) = self.get("Multigrid") {
                recommendations.push((algo, "Multigrid scales well to large problems".to_string()));
            }
            if let Some(algo) = self.get("LowStorageRK4") {
                recommendations.push((algo, "Memory-efficient time integration".to_string()));
            }
        }

        recommendations
    }

    /// Generate complexity analysis report
    pub fn generate_report(&self) -> String {
        let mut report = String::new();
        report.push_str("# CFD Algorithm Complexity Analysis Report\n\n");

        report.push_str("## Summary\n\n");
        let _ = writeln!(
            &mut report,
            "Registry contains {} algorithms\n",
            self.algorithms.len()
        );

        report.push_str("## Algorithms by Time Complexity\n\n");

        let complexities = [
            "O(1)",
            "O(log N)",
            "O(N)",
            "O(N log N)",
            "O(N²)",
            "O(N^{3/2})",
            "O(N³)",
        ];

        for complexity in &complexities {
            let algorithms: Vec<_> = self
                .algorithms
                .values()
                .filter(|algo| algo.time_complexity == *complexity)
                .collect();

            if !algorithms.is_empty() {
                let _ = writeln!(&mut report, "### {complexity}\n");
                for algo in algorithms {
                    let _ = writeln!(&mut report, "- **{}**: {}", algo.name, algo.memory_pattern);
                }
                report.push('\n');
            }
        }

        report.push_str("## Performance Recommendations\n\n");

        report.push_str("### For Small Problems (N < 10⁴)\n");
        for (algo, reason) in self.recommend_for_size(1000) {
            let _ = writeln!(&mut report, "- {}: {}", algo.name, reason);
        }

        report.push_str("\n### For Medium Problems (10⁴ ≤ N < 10⁷)\n");
        for (algo, reason) in self.recommend_for_size(100_000) {
            let _ = writeln!(&mut report, "- {}: {}", algo.name, reason);
        }

        report.push_str("\n### For Large Problems (N ≥ 10⁷)\n");
        for (algo, reason) in self.recommend_for_size(100_000_000) {
            let _ = writeln!(&mut report, "- {}: {}", algo.name, reason);
        }

        report
    }
}

impl Default for AlgorithmComplexityRegistry {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_registry_creation() {
        let registry = AlgorithmComplexityRegistry::new();
        assert!(!registry.algorithms.is_empty());
    }

    #[test]
    fn test_algorithm_lookup() {
        let registry = AlgorithmComplexityRegistry::new();
        let cg = registry.get("ConjugateGradient");
        assert!(cg.is_some());
        assert_eq!(cg.unwrap().time_complexity, "O(N^{3/2})");
    }

    #[test]
    fn test_complexity_filtering() {
        let registry = AlgorithmComplexityRegistry::new();
        let linear_algos = registry.by_time_complexity("O(N)");
        assert!(!linear_algos.is_empty());

        // Check that all returned algorithms have O(N) complexity
        for algo in linear_algos {
            assert_eq!(algo.time_complexity, "O(N)");
        }
    }

    #[test]
    fn test_recommendations() {
        let registry = AlgorithmComplexityRegistry::new();

        let small_recs = registry.recommend_for_size(1000);
        assert!(!small_recs.is_empty());

        let large_recs = registry.recommend_for_size(100_000_000);
        assert!(!large_recs.is_empty());
    }
}
