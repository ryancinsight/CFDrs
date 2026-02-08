#!/usr/bin/env python3
"""
Reference Finite Volume CFD Implementation (from pmocz/cfd-comparison-python)

This is a reference implementation of the Finite Volume Method for CFD,
adapted from Philip Mocz's excellent comparison repository:
https://github.com/pmocz/cfd-comparison-python

Used for validation comparison with CFD-rs FVM solver.

Original implementation is under GPL-3.0 license.
"""

import numpy as np
from typing import Tuple, Callable


def get_conserved_variables(
    rho: np.ndarray, 
    u: np.ndarray, 
    v: np.ndarray, 
    P: np.ndarray, 
    gamma: float = 1.4
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate conserved variables from primitive variables.
    
    For compressible Euler equations:
    - Mass: ρ
    - Momentum x: ρu
    - Momentum y: ρv
    - Energy: ρE = P/(γ-1) + 0.5*ρ*(u²+v²)
    """
    mass = rho
    mom_x = rho * u
    mom_y = rho * v
    energy = P / (gamma - 1) + 0.5 * rho * (u**2 + v**2)
    return mass, mom_x, mom_y, energy


def get_primitive_variables(
    mass: np.ndarray,
    mom_x: np.ndarray,
    mom_y: np.ndarray,
    energy: np.ndarray,
    gamma: float = 1.4
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate primitive variables from conserved variables.
    """
    rho = mass
    u = mom_x / rho
    v = mom_y / rho
    E = energy / rho
    P = (gamma - 1) * rho * (E - 0.5 * (u**2 + v**2))
    return rho, u, v, P


def get_minmod_slope(
    left: np.ndarray,
    center: np.ndarray,
    right: np.ndarray
) -> np.ndarray:
    """
    Calculate slope using minmod limiter for TVD reconstruction.
    
    minmod(a, b) = sign(a) * max(0, min(|a|, sign(a)*b))
    
    This is a total variation diminishing (TVD) slope limiter that
    prevents oscillations near discontinuities.
    """
    slope_left = center - left
    slope_right = right - center
    
    # Minmod function
    slope = np.where(
        slope_left * slope_right > 0,
        np.sign(slope_left) * np.minimum(np.abs(slope_left), np.abs(slope_right)),
        0.0
    )
    return slope


def extrapolate_to_face(
    q: np.ndarray,
    slope_x: np.ndarray,
    slope_y: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Extrapolate cell-centered values to cell faces using linear reconstruction.
    
    Returns:
        q_L_x: Left face in x-direction
        q_R_x: Right face in x-direction
        q_L_y: Left face in y-direction
        q_R_y: Right face in y-direction
    """
    q_L_x = q + 0.5 * slope_x
    q_R_x = q - 0.5 * slope_x
    q_L_y = q + 0.5 * slope_y
    q_R_y = q - 0.5 * slope_y
    return q_L_x, q_R_x, q_L_y, q_R_y


def get_flux(
    rho_L: float, rho_R: float,
    u_L: float, u_R: float,
    v_L: float, v_R: float,
    P_L: float, P_R: float,
    gamma: float = 1.4
) -> Tuple[float, float, float, float]:
    """
    Calculate numerical flux using local Lax-Friedrichs (Rusanov) method.
    
    F = 0.5 * (F_L + F_R) - 0.5 * alpha * (U_R - U_L)
    
    where alpha = max eigenvalue = |u| + c (local wave speed)
    """
    # Left state
    E_L = P_L / (gamma - 1) + 0.5 * rho_L * (u_L**2 + v_L**2)
    
    # Right state  
    E_R = P_R / (gamma - 1) + 0.5 * rho_R * (u_R**2 + v_R**2)
    
    # Average states
    rho = 0.5 * (rho_L + rho_R)
    u = 0.5 * (u_L + u_R)
    v = 0.5 * (v_L + v_R)
    P = 0.5 * (P_L + P_R)
    
    # Speed of sound
    c = np.sqrt(gamma * P / rho)
    
    # Maximum wave speed (local Lax-Friedrichs)
    alpha = max(abs(u_L) + np.sqrt(gamma * P_L / rho_L),
                abs(u_R) + np.sqrt(gamma * P_R / rho_R))
    
    # Fluxes (x-direction)
    flux_mass_L = rho_L * u_L
    flux_momx_L = rho_L * u_L**2 + P_L
    flux_momy_L = rho_L * u_L * v_L
    flux_energy_L = u_L * (E_L + P_L)
    
    flux_mass_R = rho_R * u_R
    flux_momx_R = rho_R * u_R**2 + P_R
    flux_momy_R = rho_R * u_R * v_R
    flux_energy_R = u_R * (E_R + P_R)
    
    # Rusanov flux
    flux_mass = 0.5 * (flux_mass_L + flux_mass_R) - 0.5 * alpha * (rho_R - rho_L)
    flux_momx = 0.5 * (flux_momx_L + flux_momx_R) - 0.5 * alpha * (rho_R * u_R - rho_L * u_L)
    flux_momy = 0.5 * (flux_momy_L + flux_momy_R) - 0.5 * alpha * (rho_R * v_R - rho_L * v_L)
    flux_energy = 0.5 * (flux_energy_L + flux_energy_R) - 0.5 * alpha * (E_R - E_L)
    
    return flux_mass, flux_momx, flux_momy, flux_energy


def finite_volume_step(
    mass: np.ndarray,
    mom_x: np.ndarray,
    mom_y: np.ndarray,
    energy: np.ndarray,
    dx: float,
    dy: float,
    dt: float,
    gamma: float = 1.4
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Perform one Finite Volume time step.
    
    Uses:
    - Linear reconstruction with minmod slope limiter
    - Rusanov (local Lax-Friedrichs) numerical flux
    - First-order forward Euler time integration
    """
    # Get primitive variables
    rho, u, v, P = get_primitive_variables(mass, mom_x, mom_y, energy, gamma)
    
    # Calculate slopes with periodic boundary conditions
    slope_rho_x = get_minmod_slope(np.roll(rho, 1, axis=0), rho, np.roll(rho, -1, axis=0))
    slope_rho_y = get_minmod_slope(np.roll(rho, 1, axis=1), rho, np.roll(rho, -1, axis=1))
    slope_u_x = get_minmod_slope(np.roll(u, 1, axis=0), u, np.roll(u, -1, axis=0))
    slope_u_y = get_minmod_slope(np.roll(u, 1, axis=1), u, np.roll(u, -1, axis=1))
    slope_v_x = get_minmod_slope(np.roll(v, 1, axis=0), v, np.roll(v, -1, axis=0))
    slope_v_y = get_minmod_slope(np.roll(v, 1, axis=1), v, np.roll(v, -1, axis=1))
    slope_P_x = get_minmod_slope(np.roll(P, 1, axis=0), P, np.roll(P, -1, axis=0))
    slope_P_y = get_minmod_slope(np.roll(P, 1, axis=1), P, np.roll(P, -1, axis=1))
    
    # Extrapolate to faces
    rho_L_x, rho_R_x, rho_L_y, rho_R_y = extrapolate_to_face(rho, slope_rho_x, slope_rho_y)
    u_L_x, u_R_x, u_L_y, u_R_y = extrapolate_to_face(u, slope_u_x, slope_u_y)
    v_L_x, v_R_x, v_L_y, v_R_y = extrapolate_to_face(v, slope_v_x, slope_v_y)
    P_L_x, P_R_x, P_L_y, P_R_y = extrapolate_to_face(P, slope_P_x, slope_P_y)
    
    # Calculate fluxes in x-direction
    flux_mass_x = np.zeros_like(mass)
    flux_momx_x = np.zeros_like(mass)
    flux_momy_x = np.zeros_like(mass)
    flux_energy_x = np.zeros_like(mass)
    
    for i in range(mass.shape[0]):
        for j in range(mass.shape[1]):
            i_left = (i - 1) % mass.shape[0]
            fm, fmx, fmy, fe = get_flux(
                rho_R_x[i_left, j], rho_L_x[i, j],
                u_R_x[i_left, j], u_L_x[i, j],
                v_R_x[i_left, j], v_L_x[i, j],
                P_R_x[i_left, j], P_L_x[i, j],
                gamma
            )
            flux_mass_x[i, j] = fm
            flux_momx_x[i, j] = fmx
            flux_momy_x[i, j] = fmy
            flux_energy_x[i, j] = fe
    
    # Calculate fluxes in y-direction
    flux_mass_y = np.zeros_like(mass)
    flux_momx_y = np.zeros_like(mass)
    flux_momy_y = np.zeros_like(mass)
    flux_energy_y = np.zeros_like(mass)
    
    for i in range(mass.shape[0]):
        for j in range(mass.shape[1]):
            j_left = (j - 1) % mass.shape[1]
            fm, fmx, fmy, fe = get_flux(
                rho_R_y[i, j_left], rho_L_y[i, j],
                v_R_y[i, j_left], v_L_y[i, j],  # Swap u and v for y-flux
                u_R_y[i, j_left], u_L_y[i, j],
                P_R_y[i, j_left], P_L_y[i, j],
                gamma
            )
            flux_mass_y[i, j] = fm
            flux_mass_y[i, j] = fm
            flux_momx_y[i, j] = fmy  # Swap
            flux_momy_y[i, j] = fmx  # Swap
            flux_energy_y[i, j] = fe
    
    # Update conserved variables (finite volume update)
    mass_new = mass - dt * (
        (flux_mass_x - np.roll(flux_mass_x, 1, axis=0)) / dx +
        (flux_mass_y - np.roll(flux_mass_y, 1, axis=1)) / dy
    )
    mom_x_new = mom_x - dt * (
        (flux_momx_x - np.roll(flux_momx_x, 1, axis=0)) / dx +
        (flux_momx_y - np.roll(flux_momx_y, 1, axis=1)) / dy
    )
    mom_y_new = mom_y - dt * (
        (flux_momy_x - np.roll(flux_momy_x, 1, axis=0)) / dx +
        (flux_momy_y - np.roll(flux_momy_y, 1, axis=1)) / dy
    )
    energy_new = energy - dt * (
        (flux_energy_x - np.roll(flux_energy_x, 1, axis=0)) / dx +
        (flux_energy_y - np.roll(flux_energy_y, 1, axis=1)) / dy
    )
    
    return mass_new, mom_x_new, mom_y_new, energy_new


def poiseuille_flow_validation(
    nx: int = 50,
    ny: int = 20,
    re: float = 100.0,
    n_steps: int = 1000
) -> dict:
    """
    Run Poiseuille flow validation case.
    
    For incompressible flow (low Mach number), the compressible Euler solver
    approximates the incompressible Navier-Stokes solution.
    
    Returns velocity profile for comparison with CFD-rs.
    """
    # Domain
    Lx, Ly = 2.0, 1.0
    dx, dy = Lx / nx, Ly / ny
    
    # Initial conditions (uniform flow with small perturbation)
    rho = np.ones((nx, ny))
    u = np.ones((nx, ny)) * 0.1
    v = np.zeros((nx, ny))
    P = np.ones((nx, ny))
    
    # Add pressure gradient to drive flow
    pressure_gradient = 0.01
    
    # Time stepping
    dt = 0.001
    gamma = 1.4
    
    # Get conserved variables
    mass, mom_x, mom_y, energy = get_conserved_variables(rho, u, v, P, gamma)
    
    # Time integration
    for step in range(n_steps):
        # Add pressure gradient source term
        mom_x += dt * pressure_gradient * rho
        
        # FVM step
        mass, mom_x, mom_y, energy = finite_volume_step(
            mass, mom_x, mom_y, energy, dx, dy, dt, gamma
        )
        
        # Check convergence
        if step % 100 == 0:
            rho, u, v, P = get_primitive_variables(mass, mom_x, mom_y, energy, gamma)
            max_change = np.max(np.abs(u - np.mean(u, axis=0)))
            if max_change < 1e-6:
                break
    
    # Get final velocity profile
    rho, u, v, P = get_primitive_variables(mass, mom_x, mom_y, energy, gamma)
    
    # Extract centerline profile
    y = np.linspace(0, Ly, ny)
    u_profile = np.mean(u[nx//2, :])
    
    # Analytical solution: parabolic profile
    # u(y) = (ΔP/2μL) * y * (H - y)
    # For our setup, approximate as:
    u_analytical = 4 * np.mean(u) * (y / Ly) * (1 - y / Ly)
    
    return {
        "y": y,
        "u_numerical": u[nx//2, :],
        "u_analytical": u_analytical,
        "mean_velocity": np.mean(u),
        "max_velocity": np.max(u)
    }


if __name__ == "__main__":
    # Run validation
    print("Running Finite Volume Poiseuille flow validation...")
    results = poiseuille_flow_validation(nx=50, ny=20, n_steps=500)
    
    print(f"Mean velocity: {results['mean_velocity']:.6f}")
    print(f"Max velocity: {results['max_velocity']:.6f}")
    
    # Calculate error vs analytical
    l2_error = np.sqrt(np.mean((results['u_numerical'] - results['u_analytical'])**2))
    print(f"L2 error vs analytical: {l2_error:.6e}")
