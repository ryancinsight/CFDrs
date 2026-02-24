//! Infrastructure layer — I/O adapters, storage backends, and cross-cutting
//! technical concerns that have no domain semantics.
//!
//! This layer may depend on `domain/` types but not on `application/` services.
//!
//! ## Sub-layers
//!
//! | Module | Responsibility |
//! |--------|---------------|
//! | [`io`] | STL, VTK, OpenFOAM, and Scheme mesh I/O adapters |
//! | [`storage`] | `VertexPool`, `FaceStore`, `EdgeStore`, `SlotPool`, attribute store |
//! | [`permission`] | `GhostToken`, `GhostCell`, `PermissionedArena` (branded lifetime safety) |

pub mod io;
pub mod permission;
pub mod spatial;
pub mod storage;
