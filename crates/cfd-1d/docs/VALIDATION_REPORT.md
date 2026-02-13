# cfd-1d Validation Report

## Models
- Hagen–Poiseuille (laminar circular pipes): R = 128 μ L / (π D^4); validity Re < 2000, L/D > 10.
- Darcy–Weisbach (turbulent pipes): ΔP = f (L/D) (ρ V^2/2), resistance R = (f ρ L)/(2 A D); friction factor via Colebrook–White with Haaland initial guess; laminar limit f=64/Re.

## Tests
- Analytical resistance validation for Hagen–Poiseuille.
- Laminar limit validation for Darcy–Weisbach friction factor → resistance mapping.
- Applicability via Reynolds ranges.
- Two-way branch analytical and blood rheology validation.
- Three-way branch (trifurcation) analytical split validation:
	- symmetric split target: $Q_1 = Q_2 = Q_3 = Q_0/3$
	- pressure consistency across equal daughter branches
	- Murray extension: $D_0^3 = D_1^3 + D_2^3 + D_3^3$
- Three-way branch blood-flow physiological validation:
	- daughter shear rates constrained to physiological range
	- apparent viscosities constrained to blood-reasonable bounds

## Cross-Package Comparison Workflow
- Python trifurcation cross-validation script: `validation/validate_trifurcation.py`.
- Compares `pycfdrs` trifurcation outputs against:
	- analytical Poiseuille resistance network equations,
	- Murray-law geometry constraints,
	- optional SciPy-based linear resistance solve.
- Produces a machine-readable JSON summary in repository root for traceability.

## Units and Dimensional Analysis
- Resistance units Pa·s/m^3, consistent with ΔP/Q.
- Diameter, length in meters; viscosity Pa·s; density kg/m^3; area m^2.

## Assumptions
- Incompressible Newtonian fluid; fully developed, straight pipes; appropriate roughness ranges.

## References
- Hagen (1839), Poiseuille (1840), Colebrook (1939), Haaland (1983), White (2006).
- Murray (1926), Zamir (1988), Fung (1993), Huo & Kassab (2012).
