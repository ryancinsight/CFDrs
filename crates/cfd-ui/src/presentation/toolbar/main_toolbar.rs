//! Main toolbar — file, edit, mesh, simulation action buttons.

/// Toolbar actions that can be triggered from the main toolbar.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ToolbarAction {
    // File
    NewProject,
    OpenProject,
    SaveProject,
    ImportStl,
    ExportStl,
    ExportOpenFoam,
    ExportDrawing,

    // Edit
    Undo,
    Redo,

    // Mesh
    CreateCube,
    CreateCylinder,
    CreateSphere,
    CreateCone,
    CreateTorus,
    CreatePipe,
    CsgUnion,
    CsgIntersection,
    CsgDifference,

    // View
    ShadingFlat,
    ShadingSmooth,
    ShadingWireframe,
    FitView,

    // Simulation
    SetupSimulation,
    RunSimulation,
    StopSimulation,
}
