# cfd-1d Validation Report

## Models
- Hagen–Poiseuille (laminar circular pipes): R = 128 μ L / (π D^4); validity Re < 2000, L/D > 10.
- Darcy–Weisbach (turbulent pipes): ΔP = f (L/D) (ρ V^2/2), resistance R = (f ρ L)/(2 A D); friction factor via Colebrook–White with Haaland initial guess; laminar limit f=64/Re.

## Tests
- Analytical resistance validation for Hagen–Poiseuille.
- Laminar limit validation for Darcy–Weisbach friction factor → resistance mapping.
- Applicability via Reynolds ranges.

## Units and Dimensional Analysis
- Resistance units Pa·s/m^3, consistent with ΔP/Q.
- Diameter, length in meters; viscosity Pa·s; density kg/m^3; area m^2.

## Assumptions
- Incompressible Newtonian fluid; fully developed, straight pipes; appropriate roughness ranges.

## References
- Hagen (1839), Poiseuille (1840), Colebrook (1939), Haaland (1983), White (2006).
