# Example: milestone12_report

**Source**: `crates/cfd-optim/examples/milestone12_report.rs`  
**Crate**: `cfd-optim`

## Overview

Aggregate runner over the shared Milestone 12 stage orchestration API. Accepts CLI arguments to select which optimization stages to report.

## Usage

```bash
# Report all stages
cargo run -p cfd-optim --example milestone12_report -- --all

# Report specific stages
cargo run -p cfd-optim --example milestone12_report -- --ga --option2
```

## Part Reference

Part VIII — Optimization
