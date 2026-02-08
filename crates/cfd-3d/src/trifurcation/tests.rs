#[cfg(test)]
mod tests {
    use super::super::solver::{TrifurcationSolver3D, TrifurcationConfig3D};
    use super::super::geometry::TrifurcationGeometry3D;
    use cfd_core::physics::fluid::blood::CassonBlood;
    use cfd_core::physics::fluid::traits::Fluid;
    use num_traits::Float;
    use crate::bifurcation::ConicalTransition;

    #[test]
    fn test_trifurcation_mass_conservation() {
        let l_transition = 0.002;
        let geometry = TrifurcationGeometry3D {
            d_parent: 0.01,
            l_parent: 0.05,
            d_daughters: [0.01, 0.01, 0.01],
            l_daughters: [0.02, 0.02, 0.02],
            l_transition,
            transition: ConicalTransition::SmoothCone { length: l_transition },
            branching_angles: [0.5, 0.0, -0.5], // radians (~28 degrees)
        };

        let config = TrifurcationConfig3D {
            inlet_flow_rate: 1.0e-7, // m^3/s (lower flow rate for blood)
            inlet_pressure: 100.0, // Pa
            outlet_pressures: [0.0, 0.0, 0.0],
        };

        let solver = TrifurcationSolver3D::new(geometry, config);
        
        // Use Casson blood model
        let fluid = CassonBlood::normal_blood();

        let result = solver.solve(&fluid);
        
        match result {
            Ok(solution) => {
                println!("Flow rates: {:?}", solution.flow_rates);
                println!("Mass conservation error: {}", solution.mass_conservation_error);
                
                // Check mass conservation
                // In FEM, perfect conservation is hard without very fine mesh, but should be reasonable
                assert!(solution.mass_conservation_error < 1.0e-7, "Mass conservation error too high: {}", solution.mass_conservation_error);
                
                // Check symmetry (d1 and d3 should be similar)
                let diff = Float::abs(solution.flow_rates[1] - solution.flow_rates[3]);
                assert!(diff < 1.0e-8, "Asymmetry detected: {}", diff);
            },
            Err(e) => {
                panic!("Solver failed: {}", e);
            }
        }
    }
}
