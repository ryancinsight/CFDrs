# cfd-validation

Verification and validation framework for the
[CFDrs](https://github.com/ryancinsight/CFDrs) simulation suite.

## What it provides

- **`analytical`, `analytical_benchmarks`, `solutions`** — closed-form reference
  solutions used as oracles.
- **`manufactured`** — the Method of Manufactured Solutions, for code verification
  where no analytical solution exists.
- **`convergence`** — grid and timestep refinement studies measuring observed order
  of accuracy against theoretical order.
- **`error_metrics`, `numerical`, `matrix`** — norm computation and numerical
  comparison utilities.
- **`conservation`** — mass, momentum, and energy drift checks over long integrations.
- **`literature`** — benchmarks drawn from published CFD reference cases.
- **`benchmarking`, `algorithm_complexity`, `time_integration`, `adaptive_mesh`,
  `edge_case_testing`** — performance and robustness suites.
- **`reporting`** — Markdown report generation from validation runs.

The evidence ladder is code verification (analytical/MMS), then solution verification
(discretization-error estimates), then validation against published benchmarks. A
claim is reported at the tier its evidence actually supports.

## Documentation

- [CFDrs book](https://ryancinsight.github.io/CFDrs/)
- [Workspace README](https://github.com/ryancinsight/CFDrs#readme)

## License

MIT OR Apache-2.0
