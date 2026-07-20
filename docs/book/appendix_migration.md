# Migration Notes: ndarray/nalgebra/burn → leto/hephaestus/coeus

The migration direction is:

- dense/sparse arrays and linear algebra to `leto` / `hephaestus`
- runtime parallelism to `moirai`
- SIMD surfaces to `hermes`
- ML/autodiff surfaces to `coeus`

Examples in each part are grouped to make backend replacement and validation
work incremental and testable.
