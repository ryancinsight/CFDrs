use super::manager::{TurbulenceBoundaryCondition, TurbulenceBoundaryManager};
use crate::physics::turbulence::wall_functions::{WallFunction, WallTreatment};

#[test]
fn test_wall_distance_calculation() {
    let manager = TurbulenceBoundaryManager::<f64>::new(5, 5, 0.1, 0.1);

    assert!(manager.wall_distances[0] > 0.0);
    assert!(manager.wall_distances[0] < 0.1);

    let center_idx = 2 * 5 + 2;
    assert!(manager.wall_distances[center_idx] > manager.wall_distances[0]);
}

#[test]
fn test_k_epsilon_wall_boundaries() {
    let manager = TurbulenceBoundaryManager::<f64>::new(5, 5, 0.1, 0.1);
    let boundaries = vec![(
        "south".to_string(),
        TurbulenceBoundaryCondition::Wall {
            wall_treatment: WallTreatment::new(WallFunction::Standard),
        },
    )];

    let mut k = vec![1.0; 25];
    let mut epsilon = vec![1.0; 25];

    manager.apply_k_epsilon_boundaries(&mut k, &mut epsilon, &boundaries);

    for i in 0..5 {
        assert_eq!(k[i], 0.0);
        assert!(epsilon[i] > 0.0 && epsilon[i] < 1e-6);
    }
}

#[test]
fn test_k_omega_wall_boundaries() {
    let manager = TurbulenceBoundaryManager::<f64>::new(5, 5, 0.1, 0.1);
    let boundaries = vec![(
        "south".to_string(),
        TurbulenceBoundaryCondition::Wall {
            wall_treatment: WallTreatment::new(WallFunction::Standard),
        },
    )];

    let mut k = vec![1.0; 25];
    let mut omega = vec![1.0; 25];

    manager.apply_k_omega_sst_boundaries(&mut k, &mut omega, &boundaries);

    for i in 0..5 {
        assert_eq!(k[i], 0.0);
        assert!(omega[i] > 1.0);
    }
}

#[test]
fn test_sa_wall_boundaries() {
    let manager = TurbulenceBoundaryManager::<f64>::new(5, 5, 0.1, 0.1);
    let boundaries = vec![(
        "south".to_string(),
        TurbulenceBoundaryCondition::Wall {
            wall_treatment: WallTreatment::new(WallFunction::Standard),
        },
    )];

    let mut nu_tilde = vec![1.0; 25];

    manager.apply_spalart_allmaras_boundaries(&mut nu_tilde, &boundaries);

    for i in 0..5 {
        assert_eq!(nu_tilde[i], 0.0);
    }
}

#[test]
fn test_inlet_boundaries() {
    let manager = TurbulenceBoundaryManager::<f64>::new(5, 5, 0.1, 0.1);
    let boundaries = vec![(
        "west".to_string(),
        TurbulenceBoundaryCondition::Inlet {
            turbulence_intensity: 0.05,
            turbulence_length_scale: 0.01,
            reference_velocity: 1.0,
        },
    )];

    let mut k = vec![0.0; 25];
    let mut epsilon = vec![0.0; 25];

    manager.apply_k_epsilon_boundaries(&mut k, &mut epsilon, &boundaries);

    for j in 0..5 {
        let idx = j * 5;
        assert!(k[idx] > 0.0);
        assert!(epsilon[idx] > 0.0);
    }
}
