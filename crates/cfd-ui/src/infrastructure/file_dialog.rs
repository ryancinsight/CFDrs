//! File dialog wrappers using the `rfd` crate.

use std::path::PathBuf;

/// File type filters for the dialog.
pub enum FileFilter {
    Stl,
    Obj,
    Ply,
    Glb,
    Dxf,
    Svg,
    Project,
    OpenFoam,
    /// Combined filter for all supported mesh import formats.
    AllMesh,
}

impl FileFilter {
    fn label(&self) -> &str {
        match self {
            Self::Stl => "STL Files",
            Self::Obj => "OBJ Files",
            Self::Ply => "PLY Files",
            Self::Glb => "glTF Binary Files",
            Self::Dxf => "DXF Files",
            Self::Svg => "SVG Files",
            Self::Project => "CFD Project Files",
            Self::OpenFoam => "OpenFOAM Directory",
            Self::AllMesh => "All Mesh Files",
        }
    }

    fn extensions(&self) -> &[&str] {
        match self {
            Self::Stl => &["stl"],
            Self::Obj => &["obj"],
            Self::Ply => &["ply"],
            Self::Glb => &["glb", "gltf"],
            Self::Dxf => &["dxf"],
            Self::Svg => &["svg"],
            Self::Project => &["cfdproj", "json"],
            Self::OpenFoam => &[""],
            Self::AllMesh => &["stl", "obj", "ply"],
        }
    }
}

/// Show an open-file dialog with the given filter. Returns the selected path.
pub fn open_file(filter: FileFilter) -> Option<PathBuf> {
    rfd::FileDialog::new()
        .add_filter(filter.label(), filter.extensions())
        .pick_file()
}

/// Show a save-file dialog with the given filter. Returns the selected path.
pub fn save_file(filter: FileFilter) -> Option<PathBuf> {
    rfd::FileDialog::new()
        .add_filter(filter.label(), filter.extensions())
        .save_file()
}

/// Show a folder picker dialog. Returns the selected directory path.
pub fn pick_folder() -> Option<PathBuf> {
    rfd::FileDialog::new().pick_folder()
}
