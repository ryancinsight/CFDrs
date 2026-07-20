# CFDrs Architecture and Problem Setup

CFDrs decomposes simulation concerns into focused crates (`cfd-core`,
`cfd-math`, `cfd-1d/2d/3d`, validation, and IO), then composes them through
workspace-level contracts.

The examples in this part show baseline setup, scenario wiring, and stable
entry points for expanding to larger validation workflows.
