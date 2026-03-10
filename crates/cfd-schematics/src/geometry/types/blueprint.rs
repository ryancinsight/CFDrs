// This file is intentionally empty.
//
// The `ChannelSystem::to_blueprint` conversion it previously contained has
// been consolidated into `NetworkBlueprint` directly.  All interchange,
// rendering, and JSON serialisation methods now live on
// `crate::domain::model::NetworkBlueprint`.
//
// This file is NOT part of the module tree (no `mod blueprint;` in
// `geometry/types/mod.rs`) and exists only as a historical artefact.
