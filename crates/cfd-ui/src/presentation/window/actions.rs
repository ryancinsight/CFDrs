//! Keyboard action definitions and handlers for the gpui window.

use gpui::{actions, App, Context, KeyBinding, Window};

// Define gpui action types. Each becomes a zero-sized struct that can be
// dispatched via keyboard shortcuts or toolbar clicks.
actions!(
    cfd_ui,
    [
        Undo,
        Redo,
        Delete,
        SelectAll,
        FitView,
        ToggleWireframe,
        NewProject,
        OpenProject,
        SaveProject,
        ImportStl,
        ImportMesh,
        ExportStl,
        ExportObj,
        ExportPly,
        ExportGlb,
        ExportDxf,
        ExportDrawing,
        ExportOpenFoam,
        CreateCube,
        CreateCylinder,
        CreateSphere,
        CreateCone,
        CreateTorus,
        CreatePipe,
        CreateEllipsoid,
        CreateCapsule,
        CreateFrustum,
        CreateGeodesicSphere,
        CreateRoundedCube,
        CreateElbow,
        CsgUnion,
        CsgIntersection,
        CsgDifference,
        RunSimulation,
        StopSimulation,
    ]
);

/// Register global keyboard shortcuts.
pub fn register_keybindings(app: &mut App) {
    app.bind_keys([
        KeyBinding::new("ctrl-z", Undo, None),
        KeyBinding::new("ctrl-shift-z", Redo, None),
        KeyBinding::new("delete", Delete, None),
        KeyBinding::new("ctrl-a", SelectAll, None),
        KeyBinding::new("f", FitView, None),
        KeyBinding::new("z", ToggleWireframe, None),
        KeyBinding::new("ctrl-n", NewProject, None),
        KeyBinding::new("ctrl-o", OpenProject, None),
        KeyBinding::new("ctrl-s", SaveProject, None),
    ]);
}

use super::workspace::Workspace;
use crate::application::mesh_ops::csg::{CsgBooleanCommand, CsgOp};
use crate::application::mesh_ops::import::{ImportMeshCommand, ImportStlCommand};
use crate::application::mesh_ops::primitives::{CreatePrimitiveCommand, PrimitiveSpec};
use crate::domain::scene::graph::SceneEntity;
use crate::infrastructure::gpu::pipeline::ShadingMode;

impl Workspace {
    pub(crate) fn handle_undo(
        &mut self,
        _action: &Undo,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        if let Err(e) = self.state.history.undo(&mut self.state.document) {
            self.state.console.error(format!("Undo failed: {e}"));
        }
        self.refresh_viewport();
        cx.notify();
    }

    pub(crate) fn handle_redo(
        &mut self,
        _action: &Redo,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        if let Err(e) = self.state.history.redo(&mut self.state.document) {
            self.state.console.error(format!("Redo failed: {e}"));
        }
        self.refresh_viewport();
        cx.notify();
    }

    pub(crate) fn handle_fit_view(
        &mut self,
        _action: &FitView,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        let camera = self.state.document.scene.camera_mut();
        let min = nalgebra::Point3::new(-1.0, -1.0, -1.0);
        let max = nalgebra::Point3::new(1.0, 1.0, 1.0);
        camera.fit_to_bounds(&min, &max);
        self.refresh_viewport();
        cx.notify();
    }

    pub(crate) fn handle_toggle_wireframe(
        &mut self,
        _action: &ToggleWireframe,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        if let Some(renderer) = &mut self.renderer {
            let current = renderer.shading_mode();
            let next = match current {
                ShadingMode::Wireframe => ShadingMode::Smooth,
                _ => ShadingMode::Wireframe,
            };
            renderer.set_shading_mode(next);
        }
        self.refresh_viewport();
        cx.notify();
    }

    pub(crate) fn handle_new_project(
        &mut self,
        _action: &NewProject,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.state = crate::presentation::workspace::WorkspaceState::new();
        self.state.console.info("New project created.".to_owned());
        self.refresh_viewport();
        cx.notify();
    }

    // -- Primitive creation handlers ------------------------------------------

    fn create_primitive(
        &mut self,
        spec: PrimitiveSpec,
        cx: &mut Context<Self>,
    ) {
        let name = spec.type_name().to_owned();
        let cmd = CreatePrimitiveCommand::new(spec, name.clone());
        match self.state.history.execute(Box::new(cmd), &mut self.state.document) {
            Ok(()) => self.state.console.info(format!("Created {name}")),
            Err(e) => self.state.console.error(format!("Create {name} failed: {e}")),
        }
        self.refresh_viewport();
        cx.notify();
    }

    pub(crate) fn handle_create_cube(
        &mut self,
        _action: &CreateCube,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.create_primitive(
            PrimitiveSpec::Cube { width: 1.0, height: 1.0, depth: 1.0 },
            cx,
        );
    }

    pub(crate) fn handle_create_cylinder(
        &mut self,
        _action: &CreateCylinder,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.create_primitive(
            PrimitiveSpec::Cylinder { radius: 0.5, height: 1.0, segments: 32 },
            cx,
        );
    }

    pub(crate) fn handle_create_sphere(
        &mut self,
        _action: &CreateSphere,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.create_primitive(
            PrimitiveSpec::Sphere { radius: 0.5, segments: 32, stacks: 16 },
            cx,
        );
    }

    pub(crate) fn handle_create_cone(
        &mut self,
        _action: &CreateCone,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.create_primitive(
            PrimitiveSpec::Cone { radius: 0.5, height: 1.0, segments: 32 },
            cx,
        );
    }

    pub(crate) fn handle_create_torus(
        &mut self,
        _action: &CreateTorus,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.create_primitive(
            PrimitiveSpec::Torus {
                major_radius: 1.0,
                minor_radius: 0.25,
                major_segments: 32,
                minor_segments: 16,
            },
            cx,
        );
    }

    pub(crate) fn handle_create_pipe(
        &mut self,
        _action: &CreatePipe,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.create_primitive(
            PrimitiveSpec::Pipe {
                outer_radius: 0.5,
                inner_radius: 0.4,
                height: 1.0,
                segments: 32,
            },
            cx,
        );
    }

    pub(crate) fn handle_create_ellipsoid(
        &mut self,
        _action: &CreateEllipsoid,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.create_primitive(
            PrimitiveSpec::Ellipsoid {
                semi_x: 0.5,
                semi_y: 0.35,
                semi_z: 0.25,
                segments: 32,
                stacks: 16,
            },
            cx,
        );
    }

    pub(crate) fn handle_create_capsule(
        &mut self,
        _action: &CreateCapsule,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.create_primitive(
            PrimitiveSpec::Capsule {
                radius: 0.3,
                cylinder_height: 0.6,
                segments: 32,
                hemisphere_stacks: 8,
            },
            cx,
        );
    }

    pub(crate) fn handle_create_frustum(
        &mut self,
        _action: &CreateFrustum,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.create_primitive(
            PrimitiveSpec::Frustum {
                bottom_radius: 0.5,
                top_radius: 0.25,
                height: 1.0,
                segments: 32,
            },
            cx,
        );
    }

    pub(crate) fn handle_create_geodesic_sphere(
        &mut self,
        _action: &CreateGeodesicSphere,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.create_primitive(
            PrimitiveSpec::GeodesicSphere {
                radius: 0.5,
                frequency: 3,
            },
            cx,
        );
    }

    pub(crate) fn handle_create_rounded_cube(
        &mut self,
        _action: &CreateRoundedCube,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.create_primitive(
            PrimitiveSpec::RoundedCube {
                width: 1.0,
                height: 1.0,
                depth: 1.0,
                corner_radius: 0.15,
                corner_segments: 4,
            },
            cx,
        );
    }

    pub(crate) fn handle_create_elbow(
        &mut self,
        _action: &CreateElbow,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.create_primitive(
            PrimitiveSpec::Elbow {
                tube_radius: 0.15,
                bend_radius: 0.5,
                bend_angle: std::f64::consts::FRAC_PI_2,
                tube_segments: 16,
                arc_segments: 16,
            },
            cx,
        );
    }

    // -- CSG boolean handlers -------------------------------------------------

    fn csg_op(
        &mut self,
        op: CsgOp,
        label: &str,
        cx: &mut Context<Self>,
    ) {
        let selected: Vec<usize> = self.state.selection.iter().collect();
        if selected.len() != 2 {
            self.state
                .console
                .error(format!("{label} requires exactly 2 selected meshes."));
            return;
        }
        let handle_a = match self.state.document.scene.node(selected[0]) {
            Some(n) => match &n.entity {
                SceneEntity::Mesh(h) => *h,
                _ => {
                    self.state.console.error("First selection is not a mesh.".to_owned());
                    return;
                }
            },
            None => return,
        };
        let handle_b = match self.state.document.scene.node(selected[1]) {
            Some(n) => match &n.entity {
                SceneEntity::Mesh(h) => *h,
                _ => {
                    self.state.console.error("Second selection is not a mesh.".to_owned());
                    return;
                }
            },
            None => return,
        };
        let cmd = CsgBooleanCommand::new(op, handle_a, handle_b, label.to_owned());
        match self.state.history.execute(Box::new(cmd), &mut self.state.document) {
            Ok(()) => self.state.console.info(format!("{label} completed.")),
            Err(e) => self.state.console.error(format!("{label} failed: {e}")),
        }
        self.refresh_viewport();
        cx.notify();
    }

    pub(crate) fn handle_csg_union(
        &mut self,
        _action: &CsgUnion,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.csg_op(CsgOp::Union, "CSG Union", cx);
    }

    pub(crate) fn handle_csg_intersection(
        &mut self,
        _action: &CsgIntersection,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.csg_op(CsgOp::Intersection, "CSG Intersection", cx);
    }

    pub(crate) fn handle_csg_difference(
        &mut self,
        _action: &CsgDifference,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.csg_op(CsgOp::Difference, "CSG Difference", cx);
    }

    // -- Import / Export handlers ----------------------------------------------

    pub(crate) fn handle_import_stl(
        &mut self,
        _action: &ImportStl,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        use crate::infrastructure::file_dialog::{self, FileFilter};

        let Some(path) = file_dialog::open_file(FileFilter::Stl) else {
            return;
        };
        let cmd = ImportStlCommand::new(path);
        match self.state.history.execute(Box::new(cmd), &mut self.state.document) {
            Ok(()) => self.state.console.info("STL imported.".to_owned()),
            Err(e) => self.state.console.error(format!("STL import failed: {e}")),
        }
        self.refresh_viewport();
        cx.notify();
    }

    pub(crate) fn handle_export_stl(
        &mut self,
        _action: &ExportStl,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        use crate::infrastructure::file_dialog::{self, FileFilter};

        // Export the selected mesh, or the first visible mesh.
        let handle = self.selected_mesh_handle().or_else(|| self.first_visible_mesh_handle());
        let Some(handle) = handle else {
            self.state.console.error("No mesh to export.".to_owned());
            return;
        };
        let Some(mesh) = self.state.document.mesh(handle) else {
            self.state.console.error("Mesh data not found.".to_owned());
            return;
        };
        let Some(path) = file_dialog::save_file(FileFilter::Stl) else {
            return;
        };
        match crate::application::export::stl::export_stl(mesh, &path, true) {
            Ok(()) => self.state.console.info(format!("Exported STL to {}", path.display())),
            Err(e) => self.state.console.error(format!("STL export failed: {e}")),
        }
        cx.notify();
    }

    pub(crate) fn handle_export_openfoam(
        &mut self,
        _action: &ExportOpenFoam,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        use crate::infrastructure::file_dialog;

        let handle = self.selected_mesh_handle().or_else(|| self.first_visible_mesh_handle());
        let Some(handle) = handle else {
            self.state.console.error("No mesh to export.".to_owned());
            return;
        };
        let Some(mesh) = self.state.document.mesh(handle) else {
            self.state.console.error("Mesh data not found.".to_owned());
            return;
        };
        let Some(dir) = file_dialog::pick_folder() else {
            return;
        };
        match crate::application::export::openfoam::export_openfoam(mesh, &dir, &[]) {
            Ok(()) => self.state.console.info(format!("Exported OpenFOAM to {}", dir.display())),
            Err(e) => self.state.console.error(format!("OpenFOAM export failed: {e}")),
        }
        cx.notify();
    }

    pub(crate) fn handle_import_mesh(
        &mut self,
        _action: &ImportMesh,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        use crate::infrastructure::file_dialog::{self, FileFilter};

        let Some(path) = file_dialog::open_file(FileFilter::AllMesh) else {
            return;
        };
        let cmd = ImportMeshCommand::new(path);
        match self.state.history.execute(Box::new(cmd), &mut self.state.document) {
            Ok(()) => self.state.console.info("Mesh imported.".to_owned()),
            Err(e) => self.state.console.error(format!("Mesh import failed: {e}")),
        }
        self.refresh_viewport();
        cx.notify();
    }

    pub(crate) fn handle_export_obj(
        &mut self,
        _action: &ExportObj,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.export_mesh_format("OBJ", crate::infrastructure::file_dialog::FileFilter::Obj, |mesh, path| {
            crate::application::export::obj::export_obj(mesh, path)
        });
        cx.notify();
    }

    pub(crate) fn handle_export_ply(
        &mut self,
        _action: &ExportPly,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.export_mesh_format("PLY", crate::infrastructure::file_dialog::FileFilter::Ply, |mesh, path| {
            crate::application::export::ply::export_ply(mesh, path)
        });
        cx.notify();
    }

    pub(crate) fn handle_export_glb(
        &mut self,
        _action: &ExportGlb,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.export_mesh_format("GLB", crate::infrastructure::file_dialog::FileFilter::Glb, |mesh, path| {
            crate::application::export::glb::export_glb(mesh, path)
        });
        cx.notify();
    }

    pub(crate) fn handle_export_dxf(
        &mut self,
        _action: &ExportDxf,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.export_mesh_format("DXF", crate::infrastructure::file_dialog::FileFilter::Dxf, |mesh, path| {
            crate::application::export::dxf::export_dxf(mesh, path)
        });
        cx.notify();
    }

    pub(crate) fn handle_export_drawing(
        &mut self,
        _action: &ExportDrawing,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        use crate::infrastructure::file_dialog::{self, FileFilter};

        let handle = self.selected_mesh_handle().or_else(|| self.first_visible_mesh_handle());
        let Some(handle) = handle else {
            self.state.console.error("No mesh to export.".to_owned());
            return;
        };
        let Some(mesh) = self.state.document.mesh(handle) else {
            self.state.console.error("Mesh data not found.".to_owned());
            return;
        };
        let Some(path) = file_dialog::save_file(FileFilter::Svg) else {
            return;
        };
        let sheet = crate::presentation::dialogs::drawing_dialog::standard_three_view(
            0,
            crate::domain::drawing::sheet::SheetSize::A3,
        );
        match crate::application::export::drawing_export::export_drawing_svg(&sheet, &[mesh], &path) {
            Ok(()) => self.state.console.info(format!("Exported SVG drawing to {}", path.display())),
            Err(e) => self.state.console.error(format!("SVG export failed: {e}")),
        }
        cx.notify();
    }

    // -- File management handlers ----------------------------------------------

    pub(crate) fn handle_open_project(
        &mut self,
        _action: &OpenProject,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.state.console.warn("Project serialization not yet implemented.".to_owned());
        cx.notify();
    }

    pub(crate) fn handle_save_project(
        &mut self,
        _action: &SaveProject,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.state.console.warn("Project serialization not yet implemented.".to_owned());
        cx.notify();
    }

    // -- Edit handlers --------------------------------------------------------

    pub(crate) fn handle_delete(
        &mut self,
        _action: &Delete,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        let selected: Vec<usize> = self.state.selection.iter().collect();
        if selected.is_empty() {
            return;
        }
        for &idx in &selected {
            // Remove mesh data if the node references one.
            if let Some(node) = self.state.document.scene.node(idx) {
                if let SceneEntity::Mesh(h) = &node.entity {
                    self.state.document.remove_mesh(*h);
                }
            }
            self.state.document.scene.remove_node(idx);
        }
        self.state.selection.clear();
        self.state
            .console
            .info(format!("Deleted {} object(s).", selected.len()));
        self.refresh_viewport();
        cx.notify();
    }

    pub(crate) fn handle_select_all(
        &mut self,
        _action: &SelectAll,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        let root = self.state.document.scene.root();
        let children: Vec<usize> = self.state.document.scene.children(root).to_vec();
        self.state.selection.clear();
        for idx in children {
            self.state.selection.add(idx);
        }
        cx.notify();
    }

    // -- Simulation handlers (stubs — solver logic stays in cfd-1d/2d/3d) -----

    pub(crate) fn handle_run_simulation(
        &mut self,
        _action: &RunSimulation,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.state.console.info(
            "Simulation setup required. Solver integration with cfd-2d/3d pending.".to_owned(),
        );
        cx.notify();
    }

    pub(crate) fn handle_stop_simulation(
        &mut self,
        _action: &StopSimulation,
        _window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        self.state.console.info("No simulation running.".to_owned());
        cx.notify();
    }

    // -- Helper methods -------------------------------------------------------

    /// Get the mesh handle of the single selected node, if it is a mesh.
    fn selected_mesh_handle(&self) -> Option<crate::domain::scene::graph::MeshHandle> {
        let idx = self.state.selection.single()?;
        let node = self.state.document.scene.node(idx)?;
        match &node.entity {
            SceneEntity::Mesh(h) => Some(*h),
            _ => None,
        }
    }

    /// Get the first visible mesh handle from the scene.
    fn first_visible_mesh_handle(&self) -> Option<crate::domain::scene::graph::MeshHandle> {
        self.state
            .document
            .scene
            .visible_meshes()
            .next()
            .map(|(_idx, h)| h)
    }

    /// Generic mesh export: resolve mesh, open save dialog, call exporter.
    fn export_mesh_format(
        &mut self,
        label: &str,
        filter: crate::infrastructure::file_dialog::FileFilter,
        exporter: impl FnOnce(&cfd_mesh::IndexedMesh<f64>, &std::path::Path) -> anyhow::Result<()>,
    ) {
        use crate::infrastructure::file_dialog;

        let handle = self.selected_mesh_handle().or_else(|| self.first_visible_mesh_handle());
        let Some(handle) = handle else {
            self.state.console.error("No mesh to export.".to_owned());
            return;
        };
        let Some(mesh) = self.state.document.mesh(handle) else {
            self.state.console.error("Mesh data not found.".to_owned());
            return;
        };
        let Some(path) = file_dialog::save_file(filter) else {
            return;
        };
        match exporter(mesh, &path) {
            Ok(()) => self.state.console.info(format!("Exported {label} to {}", path.display())),
            Err(e) => self.state.console.error(format!("{label} export failed: {e}")),
        }
    }
}
