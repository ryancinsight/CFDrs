use cfd_1d::{Network, NetworkBuilder, NetworkProblem, NetworkSolver};
use cfd_1d::components::channels::CircularChannel;
use cfd_1d::components::Component;
use cfd_core::physics::fluid::database::water_20c;
use proptest::prelude::*;

fn water() -> cfd_core::physics::fluid::ConstantPropertyFluid<f64> {
    water_20c::<f64>().unwrap()
}

// Property: Hagen-Poiseuille resistance is strictly monotonic with length
proptest! {
    #[test]
    fn prop_hagen_poiseuille_monotone_in_length(
        l1 in 1e-4_f64..1.0,
        l2 in 1e-4_f64..1.0,
    ) {
        let fluid = water();
        let d = 1e-3;
        
        // Sorting ensures l_short <= l_long
        let (l_short, l_long) = if l1 < l2 { (l1, l2) } else { (l2, l1) };
        
        let ch_short = CircularChannel::new(l_short, d, 0.0);
        let ch_long = CircularChannel::new(l_long, d, 0.0);

        let r_short = ch_short.resistance(&fluid);
        let r_long = ch_long.resistance(&fluid);

        if l_short < l_long {
            prop_assert!(r_short < r_long, "Resistance must increase with length");
        } else {
            prop_assert_eq!(r_short, r_long, "Equal lengths must have equal resistance");
        }
    }
}

// Property: Hagen-Poiseuille resistance is strictly monotonic inversely with diameter
proptest! {
    #[test]
    fn prop_hagen_poiseuille_monotone_in_diameter_inverse(
        d1 in 1e-5_f64..1e-2,
        d2 in 1e-5_f64..1e-2,
    ) {
        let fluid = water();
        let l = 0.1;
        
        let (d_small, d_large) = if d1 < d2 { (d1, d2) } else { (d2, d1) };
        
        let ch_small = CircularChannel::new(l, d_small, 0.0);
        let ch_large = CircularChannel::new(l, d_large, 0.0);

        let r_small = ch_small.resistance(&fluid);
        let r_large = ch_large.resistance(&fluid);

        if d_small < d_large {
            prop_assert!(r_small > r_large, "Resistance must decrease with larger diameter");
        } else {
            prop_assert_eq!(r_small, r_large, "Equal diameters must have equal resistance");
        }
    }
}

// Property: Resistance must be strictly positive and finite
proptest! {
    #[test]
    fn prop_resistance_positive_finite(
        l in 1e-6_f64..10.0,
        d in 1e-6_f64..0.1,
    ) {
        let fluid = water();
        let ch = CircularChannel::new(l, d, 0.0);
        let r = ch.resistance(&fluid);

        prop_assert!(r > 0.0, "Resistance must be strictly positive");
        prop_assert!(r.is_finite(), "Resistance must be finite");
    }
}

// Property: Parallel combinations are strictly less than any individual component
proptest! {
    #[test]
    fn prop_parallel_combination_less_than_either(
        d1 in 1e-4_f64..1e-2,
        d2 in 1e-4_f64..1e-2,
    ) {
        let fluid = water();
        let l = 0.1;
        
        let ch1 = CircularChannel::new(l, d1, 0.0);
        let ch2 = CircularChannel::new(l, d2, 0.0);
        
        let r1 = ch1.resistance(&fluid);
        let r2 = ch2.resistance(&fluid);
        
        let r_parallel = 1.0 / (1.0 / r1 + 1.0 / r2);
        
        prop_assert!(r_parallel < r1, "Parallel resistance ({}) must be < R1 ({})", r_parallel, r1);
        prop_assert!(r_parallel < r2, "Parallel resistance ({}) must be < R2 ({})", r_parallel, r2);
    }
}

// Property: Series combinations are strictly greater than any individual component
proptest! {
    #[test]
    fn prop_series_combination_greater_than_either(
        l1 in 1e-4_f64..1.0,
        l2 in 1e-4_f64..1.0,
    ) {
        let fluid = water();
        let d = 1e-3;
        
        let ch1 = CircularChannel::new(l1, d, 0.0);
        let ch2 = CircularChannel::new(l2, d, 0.0);
        
        let r1 = ch1.resistance(&fluid);
        let r2 = ch2.resistance(&fluid);
        
        let r_series = r1 + r2;
        
        prop_assert!(r_series > r1, "Series resistance ({}) must be > R1 ({})", r_series, r1);
        prop_assert!(r_series > r2, "Series resistance ({}) must be > R2 ({})", r_series, r2);
    }
}

// Property: KCL matches exactly (sum of incoming flow == sum of outgoing flow)
proptest! {
    #[test]
    fn prop_network_solver_mass_conservation(
        p_inlet in 1e3_f64..2e5_f64,
        p_outlet in 1e3_f64..1.9e5_f64,
    ) {
        // Only run if inlet pressure is noticeably higher than outlet
        prop_assume!(p_inlet > p_outlet + 100.0);

        let fluid = water();
        let mut builder = NetworkBuilder::new();
        
        let n_in = builder.add_inlet("inlet".into());
        let n_split = builder.add_junction("split".into());
        let n_out = builder.add_outlet("outlet".into());
        
        // Connect the geometry
        builder.connect_with_pipe(n_in, n_split, "pipe1".into());
        builder.connect_with_pipe(n_split, n_out, "pipe2".into());
        // And a second parallel pipe to test branching KCL
        builder.connect_with_pipe(n_split, n_out, "pipe3".into());
        
        let mut graph = builder.build().unwrap();
        
        // Define resistances dynamically
        let r1 = CircularChannel::new(0.01, 0.001, 0.0).resistance(&fluid);
        let r2 = CircularChannel::new(0.05, 0.0005, 0.0).resistance(&fluid);
        let r3 = CircularChannel::new(0.02, 0.0008, 0.0).resistance(&fluid);
        
        let edges = graph.edge_indices().collect::<Vec<_>>();
        
        if let Some(e) = graph.edge_weight_mut(edges[0]) { e.resistance = r1; }
        if let Some(e) = graph.edge_weight_mut(edges[1]) { e.resistance = r2; }
        if let Some(e) = graph.edge_weight_mut(edges[2]) { e.resistance = r3; }
        
        let mut network = Network::new(graph, fluid.clone());
        network.set_pressure(n_in, p_inlet);
        network.set_pressure(n_out, p_outlet);
        
        let problem = NetworkProblem::new(network);
        let solver = NetworkSolver::new();
        let solved = solver.solve_network(&problem).unwrap();
        
        // We know edges[0] is in->split, edges[1,2] is split->out
        // Flow = (P_source - P_target) / R
        let mut flow_in = 0.0;
        let mut flow_out = 0.0;
        
        for edge_idx in edges {
            let e = solved.graph.edge_weight(edge_idx).unwrap();
            let (src, tgt) = solved.graph.edge_endpoints(edge_idx).unwrap();
            let p_src = solved.pressures().get(&src).copied().unwrap_or(0.0);
            let p_tgt = solved.pressures().get(&tgt).copied().unwrap_or(0.0);
            let q = (p_src - p_tgt) / e.resistance;
            
            if tgt == n_split {
                flow_in += q;
            } else if src == n_split {
                flow_out += q;
            }
        }
        
        let diff = (flow_in - flow_out).abs();
        
        // Exact mass conservation to within high precision
        prop_assert!(diff < 1e-8, "KCL violated: In={}, Out={}, Diff={}", flow_in, flow_out, diff);
        prop_assert!(flow_in > 0.0, "Flow reversed");
    }
}
