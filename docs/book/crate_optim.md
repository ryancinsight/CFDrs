# Chapter 22 — Optimization

`cfd-optim` provides design-space exploration and Pareto-front optimization
for CFD geometries.  Key examples demonstrate the Latin-hypercube sampler
(backed by `tyche-core`) and the genetic-algorithm optimizer.

| Example | Demonstrates |
|---------|-------------|
| `cell_sep_audit` | Latin-hypercube audit of the cell-separation design space |
| `milestone12_validation` | Full Pareto-front validation for Milestone 12 geometry |

```bash
cargo run -p cfd-optim --example cell_sep_audit
cargo run -p cfd-optim --example milestone12_validation
```
