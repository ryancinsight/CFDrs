//! Scaling assessment and grading logic.

use super::metrics::{ScalingGrade, ScalingTestResult};

/// Assessment of scaling performance
#[derive(Debug, Clone)]
pub struct ScalingAssessment {
    /// Overall scaling grade (A, B, C, D, F)
    pub grade: ScalingGrade,
    /// Efficiency at target core count (64 cores)
    pub efficiency_at_64_cores: f64,
    /// Communication overhead at target core count
    pub comm_overhead_at_64_cores: f64,
    /// Load imbalance at target core count
    pub load_imbalance_at_64_cores: f64,
    /// Detailed assessment notes
    pub notes: Vec<String>,
    /// Recommendations for optimization
    pub recommendations: Vec<String>,
}

impl ScalingAssessment {
    /// Assess scaling performance from test results
    pub fn assess_from_results(results: &ScalingTestResult) -> Self {
        let mut assessment = Self {
            grade: ScalingGrade::F,
            efficiency_at_64_cores: 0.0,
            comm_overhead_at_64_cores: 0.0,
            load_imbalance_at_64_cores: 0.0,
            notes: Vec::new(),
            recommendations: Vec::new(),
        };

        // Find metrics for 64 cores
        if let Some(idx_64) = results.core_counts.iter().position(|&c| c == 64) {
            assessment.efficiency_at_64_cores = results
                .scaling_efficiency
                .get(idx_64)
                .copied()
                .unwrap_or(0.0);
            assessment.comm_overhead_at_64_cores = results
                .communication_trend
                .get(idx_64)
                .copied()
                .unwrap_or(0.0);
            assessment.load_imbalance_at_64_cores = results
                .load_imbalance_trend
                .get(idx_64)
                .copied()
                .unwrap_or(1.0);
        }

        // Assess overall grade
        assessment.grade = if assessment.efficiency_at_64_cores >= 0.8
            && assessment.comm_overhead_at_64_cores <= 10.0
            && assessment.load_imbalance_at_64_cores <= 1.2
        {
            ScalingGrade::A
        } else if assessment.efficiency_at_64_cores >= 0.7
            && assessment.comm_overhead_at_64_cores <= 15.0
            && assessment.load_imbalance_at_64_cores <= 1.3
        {
            ScalingGrade::B
        } else if assessment.efficiency_at_64_cores >= 0.6
            && assessment.comm_overhead_at_64_cores <= 20.0
        {
            ScalingGrade::C
        } else if assessment.efficiency_at_64_cores >= 0.5 {
            ScalingGrade::D
        } else {
            ScalingGrade::F
        };

        // Generate assessment notes
        assessment.generate_notes(results);
        assessment.generate_recommendations();

        assessment
    }

    fn generate_notes(&mut self, results: &ScalingTestResult) {
        if results.scaling_efficiency.iter().any(|&e| e < 0.5) {
            self.notes
                .push("Poor scaling efficiency detected (<50%)".to_string());
        }

        if results.communication_trend.iter().any(|&c| c > 20.0) {
            self.notes
                .push("High communication overhead (>20%)".to_string());
        }

        if results.load_imbalance_trend.iter().any(|&l| l > 1.5) {
            self.notes
                .push("Significant load imbalance detected".to_string());
        }

        match self.grade {
            ScalingGrade::A => {
                self.notes
                    .push("Excellent scaling performance - production ready".to_string());
            }
            ScalingGrade::B => {
                self.notes.push(
                    "Good scaling performance with minor optimization opportunities".to_string(),
                );
            }
            ScalingGrade::C => {
                self.notes
                    .push("Acceptable scaling but needs optimization".to_string());
            }
            ScalingGrade::D => {
                self.notes
                    .push("Poor scaling - significant improvements needed".to_string());
            }
            ScalingGrade::F => {
                self.notes
                    .push("Failing scaling - fundamental issues present".to_string());
            }
        }
    }

    fn generate_recommendations(&mut self) {
        if self.comm_overhead_at_64_cores > 15.0 {
            self.recommendations
                .push("Optimize MPI communication patterns".to_string());
            self.recommendations
                .push("Consider non-blocking communication".to_string());
        }

        if self.load_imbalance_at_64_cores > 1.3 {
            self.recommendations
                .push("Improve load balancing algorithms".to_string());
            self.recommendations
                .push("Implement dynamic load balancing".to_string());
        }

        if self.efficiency_at_64_cores < 0.7 {
            self.recommendations
                .push("Profile hotspots and optimize compute kernels".to_string());
            self.recommendations
                .push("Consider algorithmic improvements".to_string());
        }

        if self.grade == ScalingGrade::A {
            self.recommendations
                .push("Monitor performance in production".to_string());
            self.recommendations
                .push("Consider advanced optimizations for larger scales".to_string());
        }
    }
}
