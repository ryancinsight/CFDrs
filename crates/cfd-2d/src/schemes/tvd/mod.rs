//! # Total Variation Diminishing (TVD) Schemes and MUSCL Reconstruction
//!
//! ## Mathematical Foundation
//!
//! ### Total Variation Diminishing (TVD) Theory
//!
//! A numerical scheme is **Total Variation Diminishing (TVD)** if it satisfies:
//!
//! ```math
//! TV(u^{n+1}) \leq TV(u^n)
//! ```
//!
//! where the **total variation** is defined as:
//!
//! ```math
//! TV(u) = \sum_{i=-\infty}^{\infty} |u_{i+1} - u_i|
//! ```
//!
//! **Theorem (Harten, 1983)**: A conservative, monotone scheme is TVD.
//!
//! **Theorem (Harten, 1983)**: A TVD scheme with a convex flux limiter is monotonicity-preserving.
//!
//! ### Flux Limiters and TVD Constraints
//!
//! Flux limiters $\phi(r)$ must satisfy the TVD region (Sweby, 1984):
//!
//! ```math
//! \phi(r) \geq 0, \quad 0 \leq \frac{\phi(r)}{r} \leq 2, \quad 0 \leq \phi(r) \leq 2
//! ```
//!
//! where $r = \frac{u_{i+1} - u_i}{u_i - u_{i-1}}$ is the gradient ratio
//! (van Leer 1979 convention: forward / backward difference).
//!
//! ## MUSCL (Monotonic Upstream-centered Scheme for Conservation Laws)
//!
//! ### Second-Order MUSCL Reconstruction
//!
//! **MUSCL2** reconstructs interface values using limited slopes:
//!
//! ```math
//! u_{i+1/2}^L = u_i + \frac{1}{2} \phi(r_i) (u_i - u_{i-1})
//! u_{i+1/2}^R = u_{i+1} - \frac{1}{2} \phi(r_{i+1}) (u_{i+2} - u_{i+1})
//! ```
//!
//! where $\phi(r)$ is the flux limiter function.
//!
//! ### Third-Order MUSCL Reconstruction (MUSCL3)
//!
//! **MUSCL3** uses a 4-point stencil for higher accuracy:
//!
//! ```math
//! u_{i+1/2}^L = \frac{1}{6} \left[ -u_{i-1} + 5u_i + 2u_{i+1} \right] + \frac{1}{6} \phi(r_i) \left[ u_{i-1} - 3u_i + 2u_{i+1} \right]
//! ```
//!
//! ## Von Neumann Stability Analysis for TVD Schemes
//!
//! ### CFL Condition for TVD Schemes
//!
//! For explicit TVD schemes with MUSCL reconstruction, the CFL condition is:
//!
//! ```math
//! CFL = \max(|u|, |v|) \frac{\Delta t}{\Delta x} \leq 1
//! ```
//!
//! The factor depends on the specific limiter and reconstruction:
//!
//! - **MUSCL2 with minmod limiter**: CFL ≤ 1/2 (conservative)
//! - **MUSCL2 with van Leer limiter**: CFL ≤ 1
//! - **MUSCL3 with appropriate limiters**: CFL ≤ 1/2
//!
//! **Theorem (Harten, 1984)**: Under the CFL condition CFL ≤ 1, TVD schemes
//! with convex flux limiters are stable for scalar conservation laws.
//!
//! ## Local Truncation Error (LTE) Bounds
//!
//! - **Minmod limiter**: LTE = O(Δx²) (first-order near discontinuities)
//! - **van Leer limiter**: LTE = O(Δx²) (second-order in smooth regions)
//! - **Superbee limiter**: LTE = O(Δx²) (third-order in smooth regions)
//! - **MUSCL3**: LTE = O(Δx³) in smooth regions
//!
//! **Theorem (Shu, 1997)**: TVD schemes with monotone flux limiters are stable
//! under CFL condition CFL ≤ 1 for conservation laws.
//!
//! ## Flux Limiter Functions
//!
//! | Limiter | Formula | Range | Properties |
//! |---------|---------|-------|------------|
//! | Van Leer | $\phi(r) = \frac{r + \|r\|}{1 + \|r\|}$ | $[0, 2]$ | $C^1$ continuous |
//! | Van Albada | $\phi(r) = \frac{r^2 + r}{r^2 + 1}$ | $[0, 2/3)$ | $C^\infty$ continuous |
//! | Superbee | $\phi(r) = \max[0, \min(1, 2r), \min(2, r)]$ | $[0, 2]$ | Strong shock compression |
//! | MC | $\phi(r) = \max[0, \min(\frac{1+r}{2}, 2, 2r)]$ | $[0, 2]$ | Central-upwind hybrid |
//! | Minmod | $\phi(r) = \max[0, \min(1, r)]$ | $[0, 1]$ | Most diffusive |
//!
//! ## Convergence Theory
//!
//! **Theorem (Barth & Jespersen, 1989)**: MUSCL schemes achieve:
//! - **Second-order accuracy** in smooth regions with linear reconstruction
//! - **First-order accuracy** at extrema due to limiting
//!
//! **Monotonicity Preservation (Sweby, 1984)**: A scheme is monotonicity-preserving if:
//! $0 \leq \frac{\phi(r)}{r} \leq 2, \quad 0 \leq \phi(r) \leq 2$.
//! All implemented limiters satisfy this condition.
//!
//! ## QUICK Scheme
//!
//! **QUICK (Quadratic Upstream Interpolation for Convective Kinematics)** (Leonard, 1979):
//!
//! ```math
//! u_{i+1/2} = \frac{1}{8} \left[ 6u_i + 3u_{i+1} - u_{i-1} \right]
//! ```
//!
//! Third-order accurate in smooth regions; not TVD (may oscillate); CFL ≤ 0.8.
//!
//! ## References
//!
//! - **Harten, A. (1983)**. "High resolution schemes for hyperbolic conservation laws."
//!   *J. Comput. Phys.*, 49(3), 357-393.
//! - **Sweby, P.K. (1984)**. "High resolution schemes using flux limiters."
//!   *SIAM J. Numer. Anal.*, 21(5), 995-1011.
//! - **van Leer, B. (1979)**. "Towards the ultimate conservative difference scheme, V."
//!   *J. Comput. Phys.*, 32(1), 101-136.
//! - **Barth, T.J. & Jespersen, D.C. (1989)**. "Upwind schemes on unstructured meshes."
//!   *AIAA J.*, 27(9), 1260-1262.
//! - **Leonard, B.P. (1979)**. "A stable and accurate convective modelling procedure."
//!   *Comput. Methods Appl. Mech. Eng.*, 19(1), 59-98.
//! - **Shu, C.W. (1997)**. "ENO and WENO schemes for hyperbolic conservation laws."
//!   In *Advanced numerical approximation of nonlinear hyperbolic equations*, Springer.
//!
//! # Theorem
//! For scalar conservation laws, a TVD reconstruction with a flux limiter in the
//! Sweby admissible region does not create new extrema under a CFL-stable explicit
//! update.
//!
//! **Proof sketch**:
//! Harten's theorem states that a conservative monotone scheme is TVD, meaning
//! $TV(u^{n+1}) \le TV(u^n)$. Flux limiters $\phi(r)$ satisfying
//! $0 \le \phi(r) \le \min(2r, 2)$ and $\phi(1) = 1$ preserve that property for
//! the MUSCL family. This guarantee applies to the limited TVD reconstructions;
//! QUICK and WENO remain bounded only when their own nonlinear constraints hold.

mod muscl;
mod quick;

pub use muscl::{MUSCLOrder, MUSCLScheme};
pub use quick::QUICKScheme;

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Flux limiter for TVD schemes
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum FluxLimiter {
    /// No limiter (unlimited)
    None,
    /// Van Leer limiter
    VanLeer,
    /// Van Albada limiter
    VanAlbada,
    /// Superbee limiter
    Superbee,
    /// MC (Monotonized Central) limiter
    MC,
    /// Minmod limiter
    Minmod,
}

impl FluxLimiter {
    /// Apply flux limiter function
    pub fn apply<T: RealField + Copy + FromPrimitive + Copy>(&self, r: T) -> T {
        match self {
            FluxLimiter::None => T::one(),
            FluxLimiter::VanLeer => {
                if r > T::zero() {
                    let _two =
                        super::constants::to_realfield::<T>(super::constants::FLUX_LIMITER_TWO);
                    (r + r.abs()) / (T::one() + r.abs())
                } else {
                    T::zero()
                }
            }
            FluxLimiter::VanAlbada => {
                if r > T::zero() {
                    (r * r + r) / (r * r + T::one())
                } else {
                    T::zero()
                }
            }
            FluxLimiter::Superbee => {
                let two = super::constants::to_realfield::<T>(super::constants::FLUX_LIMITER_TWO);
                // Superbee limiter: φ(r) = max(0, min(1, 2r), min(2, r))
                T::zero().max(T::one().min(two * r).max(two.min(r)))
            }
            FluxLimiter::MC => {
                let two = super::constants::to_realfield::<T>(super::constants::FLUX_LIMITER_TWO);
                T::zero().max(((T::one() + r) / two).min(two.min(two * r)))
            }
            FluxLimiter::Minmod => {
                if r > T::zero() {
                    T::one().min(r)
                } else {
                    T::zero()
                }
            }
        }
    }
}
