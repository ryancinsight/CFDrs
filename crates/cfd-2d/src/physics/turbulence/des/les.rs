use super::super::traits::LESTurbulenceModel;
use super::DetachedEddySimulation;
use super::config::DESVariant;
use leto::geometry::Vector2;
use leto::{Array2, Storage, StorageMut};
use std::f64;

impl LESTurbulenceModel for DetachedEddySimulation {
    fn update(
        &mut self,
        velocity_u: &Array2<f64>,
        velocity_v: &Array2<f64>,
        _pressure: &Array2<f64>,
        _density: f64,
        _viscosity: f64,
        dt: f64,
        dx: f64,
        dy: f64,
    ) -> cfd_core::error::Result<()> {
        self.dx = dx;
        self.dy = dy;

        // Populate velocity buffer (SoA -> AoS adaptation)
        let nx = velocity_u.shape()[0];
        let ny = velocity_v.shape()[1];
        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;
                self.velocity_buffer[idx] = Vector2::new(velocity_u[[i, j]], velocity_v[[i, j]]);
            }
        }

        // Compute DES length scale
        self.des_length_scale = self.compute_des_length_scale(velocity_u, velocity_v, dx, dy);

        // Update nu_tilde using Spalart-Allmaras solver with DES length scale override
        // This effectively implements the DES formulation where the destruction term
        // uses min(d, C_DES * Delta) instead of d.
        self.spalart_allmaras.update_with_distance(
            self.nu_tilde.storage_mut().as_mut_slice(),
            &self.velocity_buffer,
            self.config.rans_viscosity, // Use config value as molecular viscosity
            self.des_length_scale.storage().as_slice(),
            dx,
            dy,
            dt,
        )?;

        // Update SGS/Turbulent viscosity field
        for j in 0..ny {
            for i in 0..nx {
                let nu_t = self
                    .spalart_allmaras
                    .eddy_viscosity(self.nu_tilde[[i, j]], self.config.rans_viscosity);

                // Limit SGS viscosity if needed (though SA usually stable)
                let max_visc = self.config.max_sgs_ratio * self.config.rans_viscosity;
                self.sgs_viscosity[[i, j]] = nu_t.min(max_visc);
            }
        }

        Ok(())
    }

    fn get_viscosity(&self, i: usize, j: usize) -> f64 {
        // Return total viscosity (molecular + SGS)
        self.config.rans_viscosity + self.sgs_viscosity[[i, j]]
    }

    fn get_turbulent_viscosity_field(&self) -> &Array2<f64> {
        // Return SGS viscosity field (turbulent viscosity only)
        // The total viscosity = molecular + turbulent is computed in get_viscosity()
        &self.sgs_viscosity
    }

    fn get_turbulent_kinetic_energy(&self, i: usize, j: usize) -> f64 {
        // For DES, TKE is typically obtained from the underlying RANS model
        // Since this implementation doesn't include RANS coupling, we estimate
        // TKE from the SGS viscosity using the relation: k ≈ (ν_sgs / C_k)^{2/3}
        // where C_k is a constant related to the SGS model

        let nu_sgs = self.sgs_viscosity[[i, j]];
        if nu_sgs > 0.0 {
            // Estimate k from SGS viscosity: k ≈ (ν_sgs / C_k)^{2/3}
            // Using C_k ≈ 0.1 (typical value for Smagorinsky-based models)
            let c_k = 0.1;
            (nu_sgs / c_k).powf(2.0 / 3.0)
        } else {
            0.0
        }
    }

    fn get_dissipation_rate(&self, i: usize, j: usize) -> f64 {
        // For DES, dissipation rate ε is related to TKE and length scale: ε = k^{3/2} / l
        // This provides a physically consistent estimate

        let k = self.get_turbulent_kinetic_energy(i, j);
        let l_des = self.des_length_scale[[i, j]];

        if k > 0.0 && l_des > 0.0 {
            k.powf(1.5) / l_des
        } else {
            0.0
        }
    }

    fn boundary_condition_update(
        &mut self,
        _boundary_manager: &super::super::boundary_conditions::TurbulenceBoundaryManager<f64>,
    ) -> cfd_core::error::Result<()> {
        // DES typically uses homogeneous boundary conditions
        // or recycling/rescaling methods for inflow
        Ok(())
    }

    fn get_model_name(&self) -> &str {
        match self.config.variant {
            DESVariant::DES97 => "DES97",
            DESVariant::DDES => "Delayed DES",
            DESVariant::IDDES => "Improved DDES",
        }
    }

    fn get_model_constants(&self) -> Vec<(&str, f64)> {
        let constants = vec![
            ("DES Constant", self.config.des_constant),
            ("Max SGS Ratio", self.config.max_sgs_ratio),
        ];

        // DES-specific constants only (RANS constants would come from separate model)

        constants
    }
}
