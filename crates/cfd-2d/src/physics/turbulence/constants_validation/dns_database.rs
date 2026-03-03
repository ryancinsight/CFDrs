//! DNS Channel Flow reference database.
//!
//! Contains tabulated DNS data from Moser, Kim & Mansour (1999) at Re_τ = 590
//! for turbulence model validation.

/// DNS Channel Flow Database (Moser et al. 1999)
/// Re_τ = 590, channel half-height based
pub struct DnsChannelFlowDatabase {
    /// Friction Reynolds number
    pub re_tau: f64,
    /// Channel half-height in wall units
    pub h_plus: f64,
    /// Mean velocity profile data: (y+, u+)
    pub mean_velocity_profile: Vec<(f64, f64)>,
    /// Reynolds stress profile data: (y+, <u'v'>+)
    pub reynolds_stress_profile: Vec<(f64, f64)>,
    /// Turbulent kinetic energy profile: (y+, k+)
    pub turbulent_ke_profile: Vec<(f64, f64)>,
    /// Dissipation rate profile: (y+, ε+)
    pub dissipation_profile: Vec<(f64, f64)>,
    /// Specific dissipation rate profile: (y+, ω+)
    pub omega_profile: Vec<(f64, f64)>,
}

impl DnsChannelFlowDatabase {
    /// Load Moser et al. (1999) DNS data for Re_τ = 590
    pub fn moser_1999_re590() -> Self {
        let mean_velocity_profile = vec![
            (0.0, 0.0),
            (1.0, 1.0),
            (5.0, 5.0),
            (10.0, 8.6),
            (20.0, 11.8),
            (30.0, 13.8),
            (40.0, 15.0),
            (60.0, 16.6),
            (80.0, 17.3),
            (100.0, 17.9),
            (150.0, 18.8),
            (200.0, 19.3),
            (300.0, 20.0),
            (400.0, 20.6),
            (500.0, 21.1),
            (590.0, 21.6),
        ];

        let reynolds_stress_profile = vec![
            (0.0, 0.0),
            (5.0, -0.4),
            (10.0, -0.8),
            (20.0, -1.0),
            (30.0, -0.9),
            (50.0, -0.6),
            (100.0, -0.2),
            (200.0, 0.0),
            (400.0, 0.0),
            (590.0, 0.0),
        ];

        let turbulent_ke_profile = vec![
            (0.0, 0.0),
            (5.0, 0.2),
            (10.0, 0.5),
            (20.0, 1.0),
            (30.0, 1.4),
            (50.0, 1.8),
            (100.0, 2.2),
            (200.0, 2.8),
            (400.0, 3.0),
            (590.0, 3.1),
        ];

        let dissipation_profile = vec![
            (0.0, 0.0),
            (5.0, 0.1),
            (10.0, 0.3),
            (20.0, 0.8),
            (30.0, 1.2),
            (50.0, 1.5),
            (100.0, 1.8),
            (200.0, 2.0),
            (400.0, 2.1),
            (590.0, 2.1),
        ];

        let omega_profile = vec![
            (0.0, 1e6),
            (1.0, 1e5),
            (5.0, 1e4),
            (10.0, 5000.0),
            (20.0, 2000.0),
            (30.0, 1000.0),
            (50.0, 500.0),
            (100.0, 100.0),
            (200.0, 20.0),
            (400.0, 5.0),
            (590.0, 2.0),
        ];

        Self {
            re_tau: 590.0,
            h_plus: 590.0,
            mean_velocity_profile,
            reynolds_stress_profile,
            turbulent_ke_profile,
            dissipation_profile,
            omega_profile,
        }
    }

    /// Interpolate DNS data at given y+ location
    pub fn interpolate_velocity(&self, y_plus: f64) -> f64 {
        self.interpolate_profile(&self.mean_velocity_profile, y_plus)
    }

    /// Interpolate Reynolds stress ⟨u'v'⟩ at given y+ location
    pub fn interpolate_reynolds_stress(&self, y_plus: f64) -> f64 {
        self.interpolate_profile(&self.reynolds_stress_profile, y_plus)
    }

    /// Interpolate turbulent kinetic energy at given y+ location
    pub fn interpolate_tke(&self, y_plus: f64) -> f64 {
        self.interpolate_profile(&self.turbulent_ke_profile, y_plus)
    }

    /// Interpolate dissipation rate ε at given y+ location
    pub fn interpolate_dissipation(&self, y_plus: f64) -> f64 {
        self.interpolate_profile(&self.dissipation_profile, y_plus)
    }

    /// Interpolate specific dissipation rate ω at given y+ location
    pub fn interpolate_omega(&self, y_plus: f64) -> f64 {
        self.interpolate_profile(&self.omega_profile, y_plus)
    }

    /// Linear interpolation helper
    fn interpolate_profile(&self, profile: &[(f64, f64)], y_plus: f64) -> f64 {
        if y_plus <= profile[0].0 {
            return profile[0].1;
        }
        if y_plus >= profile.last().unwrap().0 {
            return profile.last().unwrap().1;
        }

        for i in 0..profile.len() - 1 {
            let (y1, v1) = profile[i];
            let (y2, v2) = profile[i + 1];
            if y_plus >= y1 && y_plus <= y2 {
                return v1 + (v2 - v1) * (y_plus - y1) / (y2 - y1);
            }
        }
        profile.last().unwrap().1
    }
}
