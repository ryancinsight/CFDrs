#!/usr/bin/env python3
"""
External Reference Implementation: Lid-Driven Cavity Flow

Pure Python implementation using finite difference method for comparison with CFD-RS.
Based on:
- Ghia et al. (1982) "High-Re solutions for incompressible flow using Navier-Stokes equations"
- Python_CFD tutorials
- pmocz finite volume methods

Physics:
    ∂u/∂t + u·∇u = -∇p/ρ + ν∇²u    (Momentum)
    ∇·u = 0                          (Continuity)

Boundary Conditions:
    - Top wall: u = U_lid, v = 0 (moving lid)
    - Other walls: u = v = 0 (no-slip)

Validation Metrics:
    - Centerline velocities (u along vertical, v along horizontal)
    - Vortex center location
    - Minimum stream function value
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, Dict

class CavityFlowSolver:
    """
    Lid-driven cavity flow solver using projection method.
    
    Algorithm:
        1. Predict velocity (explicit advection + diffusion)
        2. Solve Poisson equation for pressure
        3. Correct velocity to satisfy continuity
        4. Repeat until convergence
    """
    
    def __init__(self, nx: int, ny: int, Re: float, dt: float = None):
        """
        Initialize cavity flow solver.
        
        Args:
            nx, ny: Grid points in x, y directions
            Re: Reynolds number (based on cavity length and lid velocity)
            dt: Time step (auto-computed if None)
        """
        self.nx = nx
        self.ny = ny
        self.Re = Re
        
        # Physical domain: [0,1] x [0,1]
        self.Lx = 1.0
        self.Ly = 1.0
        self.dx = self.Lx / (nx - 1)
        self.dy = self.Ly / (ny - 1)
        
        # Grid
        self.x = np.linspace(0, self.Lx, nx)
        self.y = np.linspace(0, self.Ly, ny)
        self.X, self.Y = np.meshgrid(self.x, self.y)
        
        # Kinematic viscosity
        self.nu = 1.0 / Re
        
        # Time step (CFL condition)
        if dt is None:
            U_lid = 1.0
            self.dt = 0.5 * min(self.dx, self.dy) / U_lid  # CFL
            self.dt = min(self.dt, 0.125 * min(self.dx, self.dy)**2 / self.nu)  # Diffusion
        else:
            self.dt = dt
        
        # Initialize fields
        self.u = np.zeros((ny, nx))  # x-velocity
        self.v = np.zeros((ny, nx))  # y-velocity
        self.p = np.zeros((ny, nx))  # pressure
        
        # Lid velocity
        self.U_lid = 1.0
        
        print(f"Cavity Flow Solver Initialized:")
        print(f"  Grid: {nx} x {ny}")
        print(f"  Re = {Re}")
        print(f"  nu = {self.nu:.6f}")
        print(f"  dx = {self.dx:.6f}, dy = {self.dy:.6f}")
        print(f"  dt = {self.dt:.6f}")
    
    def apply_boundary_conditions(self):
        """Apply no-slip and lid boundary conditions"""
        # Bottom wall: u = v = 0
        self.u[0, :] = 0.0
        self.v[0, :] = 0.0
        
        # Top wall: u = U_lid, v = 0
        self.u[-1, :] = self.U_lid
        self.v[-1, :] = 0.0
        
        # Left wall: u = v = 0
        self.u[:, 0] = 0.0
        self.v[:, 0] = 0.0
        
        # Right wall: u = v = 0
        self.u[:, -1] = 0.0
        self.v[:, -1] = 0.0
    
    def build_pressure_poisson_matrix(self) -> np.ndarray:
        """
        Build coefficient matrix for pressure Poisson equation.
        
        Equation: ∇²p = ρ/dt * ∇·u*
        
        Discretization (5-point stencil):
            (p[i+1,j] - 2p[i,j] + p[i-1,j]) / dx² + 
            (p[i,j+1] - 2p[i,j] + p[i,j-1]) / dy² = b[i,j]
        
        Returns:
            Sparse matrix A where A @ p_vec = b_vec
        """
        N = (self.nx - 2) * (self.ny - 2)  # Interior points
        A = np.zeros((N, N))
        
        idx = lambda i, j: i * (self.nx - 2) + j
        
        for i in range(1, self.ny - 1):
            for j in range(1, self.nx - 1):
                ii = idx(i - 1, j - 1)
                
                # Diagonal
                A[ii, ii] = -2.0 / self.dx**2 - 2.0 / self.dy**2
                
                # Off-diagonals
                if i > 1:  # Bottom neighbor
                    A[ii, idx(i - 2, j - 1)] = 1.0 / self.dy**2
                if i < self.ny - 2:  # Top neighbor
                    A[ii, idx(i, j - 1)] = 1.0 / self.dy**2
                if j > 1:  # Left neighbor
                    A[ii, idx(i - 1, j - 2)] = 1.0 / self.dx**2
                if j < self.nx - 2:  # Right neighbor
                    A[ii, idx(i - 1, j)] = 1.0 / self.dx**2
        
        return A
    
    def solve_pressure_poisson(self, u_star: np.ndarray, v_star: np.ndarray, 
                              max_iter: int = 1000, tol: float = 1e-6):
        """
        Solve pressure Poisson equation using iterative method.
        
        ∇²p = ρ/dt * (∂u*/∂x + ∂v*/∂y)
        
        Uses Gauss-Seidel iteration with Neumann BCs (∂p/∂n = 0 on walls).
        """
        p_new = self.p.copy()
        
        # RHS: divergence of predicted velocity
        rhs = np.zeros_like(self.p)
        rhs[1:-1, 1:-1] = (
            (u_star[1:-1, 2:] - u_star[1:-1, :-2]) / (2 * self.dx) +
            (v_star[2:, 1:-1] - v_star[:-2, 1:-1]) / (2 * self.dy)
        ) / self.dt
        
        # Iterative solver
        for iteration in range(max_iter):
            p_old = p_new.copy()
            
            # Interior points (Gauss-Seidel)
            for i in range(1, self.ny - 1):
                for j in range(1, self.nx - 1):
                    p_new[i, j] = (
                        (self.dy**2 * (p_new[i, j+1] + p_new[i, j-1]) +
                         self.dx**2 * (p_new[i+1, j] + p_new[i-1, j]) -
                         self.dx**2 * self.dy**2 * rhs[i, j])
                        / (2 * (self.dx**2 + self.dy**2))
                    )
            
            # Neumann BCs (∂p/∂n = 0)
            p_new[0, :] = p_new[1, :]      # Bottom
            p_new[-1, :] = p_new[-2, :]    # Top
            p_new[:, 0] = p_new[:, 1]      # Left
            p_new[:, -1] = p_new[:, -2]    # Right
            
            # Check convergence
            res = np.linalg.norm(p_new - p_old)
            if res < tol:
                break
        
        self.p = p_new
        return iteration + 1
    
    def step(self) -> float:
        """
        Advance one time step using fractional step method.
        
        Returns:
            Maximum velocity residual
        """
        u_old = self.u.copy()
        v_old = self.v.copy()
        
        # Step 1: Predict velocity (advection + diffusion, no pressure)
        u_star = u_old.copy()
        v_star = v_old.copy()
        
        # Interior points
        for i in range(1, self.ny - 1):
            for j in range(1, self.nx - 1):
                # u-momentum
                u_adv = (
                    u_old[i, j] * (u_old[i, j+1] - u_old[i, j-1]) / (2 * self.dx) +
                    v_old[i, j] * (u_old[i+1, j] - u_old[i-1, j]) / (2 * self.dy)
                )
                u_diff = self.nu * (
                    (u_old[i, j+1] - 2*u_old[i, j] + u_old[i, j-1]) / self.dx**2 +
                    (u_old[i+1, j] - 2*u_old[i, j] + u_old[i-1, j]) / self.dy**2
                )
                u_star[i, j] = u_old[i, j] + self.dt * (-u_adv + u_diff)
                
                # v-momentum
                v_adv = (
                    u_old[i, j] * (v_old[i, j+1] - v_old[i, j-1]) / (2 * self.dx) +
                    v_old[i, j] * (v_old[i+1, j] - v_old[i-1, j]) / (2 * self.dy)
                )
                v_diff = self.nu * (
                    (v_old[i, j+1] - 2*v_old[i, j] + v_old[i, j-1]) / self.dx**2 +
                    (v_old[i+1, j] - 2*v_old[i, j] + v_old[i-1, j]) / self.dy**2
                )
                v_star[i, j] = v_old[i, j] + self.dt * (-v_adv + v_diff)
        
        # Apply BCs to predicted velocity
        u_star[0, :] = 0.0; u_star[-1, :] = self.U_lid
        u_star[:, 0] = 0.0; u_star[:, -1] = 0.0
        v_star[0, :] = 0.0; v_star[-1, :] = 0.0
        v_star[:, 0] = 0.0; v_star[:, -1] = 0.0
        
        # Step 2: Solve pressure Poisson
        self.solve_pressure_poisson(u_star, v_star)
        
        # Step 3: Correct velocity
        for i in range(1, self.ny - 1):
            for j in range(1, self.nx - 1):
                self.u[i, j] = u_star[i, j] - self.dt * (self.p[i, j+1] - self.p[i, j-1]) / (2 * self.dx)
                self.v[i, j] = v_star[i, j] - self.dt * (self.p[i+1, j] - self.p[i-1, j]) / (2 * self.dy)
        
        # Apply BCs
        self.apply_boundary_conditions()
        
        # Compute residual
        residual = max(np.max(np.abs(self.u - u_old)), np.max(np.abs(self.v - v_old)))
        return residual
    
    def solve(self, max_steps: int = 10000, tol: float = 1e-6, print_interval: int = 100) -> Dict:
        """
        Solve until steady state.
        
        Returns:
            Dictionary with solution details
        """
        print(f"\nSolving cavity flow (Re={self.Re})...")
        
        for step in range(max_steps):
            res = self.step()
            
            if (step + 1) % print_interval == 0:
                print(f"  Step {step+1:5d}: residual = {res:.6e}")
            
            if res < tol:
                print(f"  Converged at step {step+1} (residual = {res:.6e})")
                break
        
        # Extract centerline velocities
        u_centerline_vertical = self.u[:, self.nx // 2]
        v_centerline_horizontal = self.v[self.ny // 2, :]
        
        # Find vortex center (min stream function ≈ min pressure for steady flow)
        i_vortex, j_vortex = np.unravel_index(np.argmin(self.p[1:-1, 1:-1]), self.p[1:-1, 1:-1].shape)
        x_vortex = self.x[j_vortex + 1]
        y_vortex = self.y[i_vortex + 1]
        
        return {
            "converged": res < tol,
            "steps": step + 1,
            "residual": res,
            "u_centerline": u_centerline_vertical,
            "v_centerline": v_centerline_horizontal,
            "vortex_center": (x_vortex, y_vortex),
            "u": self.u,
            "v": self.v,
            "p": self.p,
            "x": self.x,
            "y": self.y
        }
    
    def compare_with_ghia(self, solution: Dict) -> Dict:
        """
        Compare with Ghia et al. (1982) benchmark data.
        
        Ghia provides centerline velocities for Re = 100, 400, 1000, 3200, 5000, 7500, 10000
        """
        # Ghia Re=100 data (extracted from paper)
        ghia_re100_u_vertical = {
            "y": np.array([1.0000, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000,
                          0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000]),
            "u": np.array([1.0000, 0.8412, 0.7887, 0.7372, 0.6872, 0.2315, 0.0033, -0.1364,
                          -0.2058, -0.3109, -0.4255, -0.5264, -0.4091, -0.3708, -0.3307, 0.0000])
        }
        
        ghia_re100_v_horizontal = {
            "x": np.array([1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047,
                          0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000]),
            "v": np.array([0.0000, -0.0547, -0.0618, -0.0709, -0.0813, -0.1025, -0.1623, -0.2260,
                          0.0258, 0.3284, 0.3293, 0.3715, 0.3330, 0.3195, 0.3134, 0.3086, 0.0000])
        }
        
        if abs(self.Re - 100) < 1:
            # Interpolate solution to Ghia points
            u_interp = np.interp(ghia_re100_u_vertical["y"], self.y, solution["u_centerline"])
            u_error = np.mean(np.abs(u_interp - ghia_re100_u_vertical["u"]))
            
            v_interp = np.interp(ghia_re100_v_horizontal["x"], self.x, solution["v_centerline"])
            v_error = np.mean(np.abs(v_interp - ghia_re100_v_horizontal["v"]))
            
            return {
                "re": 100,
                "u_error": u_error,
                "v_error": v_error,
                "ghia_u": ghia_re100_u_vertical,
                "ghia_v": ghia_re100_v_horizontal
            }
        
        return {"re": self.Re, "u_error": None, "v_error": None}


def run_cavity_validation(Re: float = 100, nx: int = 65, ny: int = 65) -> Dict:
    """
    Run cavity flow validation and return results.
    
    Args:
        Re: Reynolds number
        nx, ny: Grid resolution
    
    Returns:
        Dictionary with solution and validation metrics
    """
    solver = CavityFlowSolver(nx, ny, Re)
    solution = solver.solve(max_steps=20000, tol=1e-6, print_interval=500)
    
    # Compare with Ghia benchmark
    validation = solver.compare_with_ghia(solution)
    
    return {
        "solution": solution,
        "validation": validation,
        "solver": solver
    }


if __name__ == "__main__":
    print("="*70)
    print("EXTERNAL REFERENCE: LID-DRIVEN CAVITY FLOW")
    print("="*70)
    
    # Run Re=100 case
    result = run_cavity_validation(Re=100, nx=65, ny=65)
    
    print(f"\nValidation against Ghia et al. (1982):")
    print(f"  Re = {result['validation']['re']}")
    if result['validation']['u_error'] is not None:
        print(f"  U centerline error: {result['validation']['u_error']:.6f}")
        print(f"  V centerline error: {result['validation']['v_error']:.6f}")
    
    print(f"\nVortex center: ({result['solution']['vortex_center'][0]:.4f}, {result['solution']['vortex_center'][1]:.4f})")
    
    # Plot
    try:
        import matplotlib.pyplot as plt
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Velocity magnitude
        U_mag = np.sqrt(result['solution']['u']**2 + result['solution']['v']**2)
        im0 = axes[0, 0].contourf(result['solver'].X, result['solver'].Y, U_mag, levels=20, cmap='viridis')
        axes[0, 0].set_title('Velocity Magnitude')
        axes[0, 0].set_xlabel('x')
        axes[0, 0].set_ylabel('y')
        plt.colorbar(im0, ax=axes[0, 0])
        
        # Pressure
        im1 = axes[0, 1].contourf(result['solver'].X, result['solver'].Y, result['solution']['p'], levels=20, cmap='RdBu_r')
        axes[0, 1].set_title('Pressure')
        axes[0, 1].set_xlabel('x')
        axes[0, 1].set_ylabel('y')
        plt.colorbar(im1, ax=axes[0, 1])
        
        # U centerline comparison
        axes[1, 0].plot(result['solution']['u_centerline'], result['solver'].y, 'b-', label='Computed', linewidth=2)
        if result['validation']['u_error'] is not None:
            axes[1, 0].plot(result['validation']['ghia_u']['u'], result['validation']['ghia_u']['y'], 
                           'ro', label='Ghia et al. (1982)', markersize=6)
        axes[1, 0].set_xlabel('u')
        axes[1, 0].set_ylabel('y')
        axes[1, 0].set_title('U Centerline (x=0.5)')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # V centerline comparison
        axes[1, 1].plot(result['solver'].x, result['solution']['v_centerline'], 'b-', label='Computed', linewidth=2)
        if result['validation']['v_error'] is not None:
            axes[1, 1].plot(result['validation']['ghia_v']['x'], result['validation']['ghia_v']['v'], 
                           'ro', label='Ghia et al. (1982)', markersize=6)
        axes[1, 1].set_xlabel('x')
        axes[1, 1].set_ylabel('v')
        axes[1, 1].set_title('V Centerline (y=0.5)')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('external_cavity_re100_validation.png', dpi=150, bbox_inches='tight')
        print(f"\nPlot saved: external_cavity_re100_validation.png")
        
    except ImportError:
        print("\nMatplotlib not available for plotting")
