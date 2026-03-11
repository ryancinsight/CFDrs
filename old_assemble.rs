
> fn assemble_fluid_mesh(meshes: Vec<IndexedMesh>) -> MeshResult<IndexedMesh> {
      let mut iter = meshes.into_iter();
      let first = iter.next().ok_or_else(|| MeshError::ChannelError {
          message: "no segment meshes to assemble".to_string(),
      })?;
      let mut accumulated = first;
      for mesh in iter {
          accumulated = csg_boolean_indexed(BooleanOp::Union, &accumulated, &mesh)?;
      }
  
      // Post-CSG repair: the iterative union can leave a small number of
      // misoriented faces (winding flips at T-junction seams) and phantom
      // disconnected triangles.  orient_outward() corrects the winding via
      // BFS flood-fill, and retain_largest_component() removes any floating
      // island faces that individually pass watertight=true but corrupt the
      // Euler characteristic of the whole mesh.
      accumulated.orient_outward();
      accumulated.retain_largest_component();
      accumulated.rebuild_edges();
  
      if accumulated.signed_volume() < 0.0 {
          accumulated.flip_faces();
      }
  
      if !accumulated.is_watertight() {
          let count = accumulated
              .edges_ref()
              .map_or(0, |e| e.boundary_edges().len());
          return Err(MeshError::NotWatertight { count });
      }
      Ok(accumulated)
  }
  
  // ΓöÇΓöÇ Venturi chain concatenated sweep (no CSG) ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ
  
  /// Build the fluid mesh for a VenturiChain topology ΓÇö NO CSG union required.
  ///
  /// ## Algorithm
  ///
  /// For a venturi with segments `[seg_0 (R_in), seg_1 (R_throat), seg_2 (R_in)]`, the
  /// assembled mesh consists of:
  ///
  /// 1. **Start cap** ΓÇö outward-facing fan for `seg_0`.
  /// 2. **Lateral surface** per segment ΓÇö quad-strip connecting adjacent rings.
  /// 3. **Annular cap** at each cross-section change ΓÇö connects the outer ring
  ///    (`R_large`) to the inner ring (`R_small`) at the same axial position.
  ///    This ring-pair forms a closed annular disk.
  /// 4. **End cap** ΓÇö outward-facing fan for the last segment.
  ///
  /// ## Why not CSG union?
  ///

