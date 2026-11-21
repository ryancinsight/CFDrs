# ADR 0004: Type junction geometry and minor loss

- Status: Proposed
- Class: [major]
- Date: 2026-11-21

Junction lengths, areas, density, and reference flow have distinct
Aequitas dimensions at the model's scalar precision. Existing serde keys
retain their SI suffixes and encode scalars. Rust callers use the renamed
quantity fields and extract the scalar only at the resistance trait's
assembly boundary.

From Bernoulli's dynamic head and `V = Q/A`, the signed loss is
`ΔP = C Q |Q|`, where `C = Kρ/(2A²)`. The standalone correlation and
blueprint network use the same coefficient calculation with their own
K-factor laws. The configured density applies to the reference derivative,
while fluid properties supply the density for a solved junction.
The derivative `2C|Q|` is even and has hydraulic resistance dimensions.

Explicit validation and coefficient evaluation reject nonfinite or
out-of-domain geometry. Keeping the existing wire keys preserves
saved model configurations. Tests verify the area ratio, both
flow signs, and exact numerical coefficients and serialized values.
