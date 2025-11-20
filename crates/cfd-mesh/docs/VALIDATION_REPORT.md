# cfd-mesh Validation Report

## Models
- Mesh topology (vertex/edge/face/cell), connectivity, boundary markers, statistics.
- Quality metrics: Knupp(2001) algebraic Jacobian q = |det J| / (||J||_F ||J^{-1}||_F), skewness CV(r=||x-C||), aspect ratio, orthogonality proxy.

## Tests
- Topology operations and statistics tested; boundary marking verified.
- Quality: unit Tet Jacobian q=1.0 (affine ref), degenerate flat q<1e-10, regular Tet skew≈0.408 (known CV).
- Planned geometry computations (areas/volumes) with analytical validation.

## Invariants
- Index validity for edges/faces/cells; consistent element type counts; zero-copy iterators for faces/vertices.
- Knupp q invariant to rigid/scaling/affine; ordered_element_vertices canonical idx order.

## Units
- Coordinates in meters; dimensional consistency for geometry when implemented.

## Assumptions
- Regular indexing; element construction APIs dictate topology correctness; sequential vertex idx canonical (v0 apex).

## References
- Knupp, P.M. (2001). Algebraic mesh quality metrics. Sandia Report SAND2001-3477.
- Shewchuk, J.R. (2002). What is a good linear finite element?
- Standard computational geometry practices for CFD meshes.
