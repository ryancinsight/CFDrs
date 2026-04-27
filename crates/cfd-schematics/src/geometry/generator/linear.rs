use crate::domain::model::{ChannelShape, ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use crate::geometry::metadata::{
    GeometryAuthoringProvenance, NodeLayoutMetadata, VenturiGeometryMetadata,
};
use crate::topology::{BlueprintTopologySpec, ParallelChannelSpec, SeriesChannelSpec};
use crate::visualizations::schematic::materialize_blueprint_layout;

fn make_outline(box_dims: (f64, f64)) -> Vec<((f64, f64), (f64, f64))> {
    let (w, h) = box_dims;
    vec![
        ((0.0, 0.0), (w, 0.0)),
        ((w, 0.0), (w, h)),
        ((w, h), (0.0, h)),
        ((0.0, h), (0.0, 0.0)),
    ]
}

fn layout_node(id: impl Into<String>, kind: NodeKind, point: (f64, f64)) -> NodeSpec {
    NodeSpec::new_at(id, kind, point).with_layout(NodeLayoutMetadata {
        x_mm: point.0,
        y_mm: point.1,
    })
}

fn serpentine_path(
    start: (f64, f64),
    end: (f64, f64),
    segments: usize,
    amplitude_mm: f64,
) -> Vec<(f64, f64)> {
    let lobes = segments.max(2);
    let dx = end.0 - start.0;
    let dy = end.1 - start.1;
    let path_len = lobes * 2;
    let mut path = Vec::with_capacity(path_len + 1);
    path.push(start);
    for idx in 1..path_len {
        let t = idx as f64 / path_len as f64;
        let x = start.0 + dx * t;
        let y = start.1
            + dy * t
            + if idx % 2 == 0 {
                amplitude_mm
            } else {
                -amplitude_mm
            };
        path.push((x, y));
    }
    path.push(end);
    path
}

fn apply_series_channel_metadata(
    mut channel: ChannelSpec,
    route: &crate::topology::ChannelRouteSpec,
    venturi: Option<&VenturiGeometryMetadata>,
) -> ChannelSpec {
    channel.therapy_zone = Some(route.therapy_zone);
    if let Some(serpentine) = &route.serpentine {
        channel.channel_shape = ChannelShape::Serpentine {
            segments: serpentine.segments,
            bend_radius_m: serpentine.bend_radius_m,
            wave_type: serpentine.wave_type,
        };
    }
    if let Some(venturi) = venturi {
        channel = channel.with_venturi_geometry(venturi.clone());
    }
    channel
}

fn series_channel_venturi(
    spec: &BlueprintTopologySpec,
    channel_id: &str,
) -> Option<VenturiGeometryMetadata> {
    spec.venturi_placements
        .iter()
        .find(|placement| placement.target_channel_id == channel_id)
        .map(|placement| VenturiGeometryMetadata {
            throat_width_m: placement.throat_geometry.throat_width_m,
            throat_height_m: placement.throat_geometry.throat_height_m,
            throat_length_m: placement.throat_geometry.throat_length_m,
            inlet_width_m: placement.throat_geometry.inlet_width_m,
            outlet_width_m: placement.throat_geometry.outlet_width_m,
            convergent_half_angle_deg: placement.throat_geometry.convergent_half_angle_deg,
            divergent_half_angle_deg: placement.throat_geometry.divergent_half_angle_deg,
            throat_position: 0.5,
        })
}

fn build_series_channel(
    from: &str,
    to: &str,
    start: (f64, f64),
    end: (f64, f64),
    channel: &SeriesChannelSpec,
    spec: &BlueprintTopologySpec,
) -> ChannelSpec {
    let mut built = ChannelSpec::new_pipe_rect(
        &channel.channel_id,
        from,
        to,
        channel.route.length_m,
        channel.route.width_m,
        channel.route.height_m,
        0.0,
        0.0,
    );
    built.path = if let Some(serpentine) = &channel.route.serpentine {
        serpentine_path(
            start,
            end,
            serpentine.segments,
            (serpentine.bend_radius_m * 1.0e3).max(channel.route.width_m * 1.5e3),
        )
    } else {
        vec![start, end]
    };
    apply_series_channel_metadata(
        built,
        &channel.route,
        series_channel_venturi(spec, &channel.channel_id).as_ref(),
    )
}

pub fn create_series_geometry_from_spec(spec: &BlueprintTopologySpec) -> NetworkBlueprint {
    let channel_count = spec.series_channels.len();
    let y_mm = spec.box_dims_mm.1 * 0.5;

    let mut nodes = Vec::with_capacity(channel_count + 1);
    let mut current_x = 10.0;
    nodes.push(layout_node("inlet", NodeKind::Inlet, (current_x, y_mm)));

    for (idx, channel) in spec.series_channels.iter().enumerate() {
        if idx < channel_count - 1 {
            current_x += channel.route.length_m * 1000.0;
            nodes.push(layout_node(
                format!("junction_{idx}"),
                NodeKind::Junction,
                (current_x, y_mm),
            ));
        } else {
            current_x += channel.route.length_m * 1000.0;
            nodes.push(layout_node("outlet", NodeKind::Outlet, (current_x, y_mm)));
        }
    }

    let mut channels = Vec::with_capacity(channel_count);
    for (idx, channel) in spec.series_channels.iter().enumerate() {
        let from = if idx == 0 {
            "inlet".to_string()
        } else {
            format!("junction_{}", idx - 1)
        };
        let to = if idx + 1 == channel_count {
            "outlet".to_string()
        } else {
            format!("junction_{idx}")
        };
        let start = nodes
            .iter()
            .find(|node| node.id.as_str() == from)
            .expect("series start node must exist")
            .point;
        let end = nodes
            .iter()
            .find(|node| node.id.as_str() == to)
            .expect("series end node must exist")
            .point;
        channels.push(build_series_channel(&from, &to, start, end, channel, spec));
    }

    let mut blueprint = NetworkBlueprint {
        name: spec.design_name.clone(),
        box_dims: spec.box_dims_mm,
        box_outline: make_outline(spec.box_dims_mm),
        nodes,
        channels,
        render_hints: None,
        topology: Some(spec.clone()),
        lineage: None,
        metadata: None,
        geometry_authored: true,
    };
    blueprint.insert_metadata(GeometryAuthoringProvenance::create_geometry());
    materialize_blueprint_layout(&mut blueprint);
    blueprint
}

fn build_parallel_channel(
    channel: &ParallelChannelSpec,
    inlet: &NodeSpec,
    outlet: &NodeSpec,
    lane_y_mm: f64,
) -> ChannelSpec {
    let mut built = ChannelSpec::new_pipe_rect(
        &channel.channel_id,
        inlet.id.as_str(),
        outlet.id.as_str(),
        channel.route.length_m,
        channel.route.width_m,
        channel.route.height_m,
        0.0,
        0.0,
    );
    let start = inlet.point;
    let end = outlet.point;
    built.path = if let Some(serpentine) = &channel.route.serpentine {
        let mut path = serpentine_path(
            (start.0 + 10.0, lane_y_mm),
            (end.0 - 10.0, lane_y_mm),
            serpentine.segments,
            (serpentine.bend_radius_m * 1.0e3).max(channel.route.width_m * 1.5e3),
        );
        path.insert(0, start);
        path.push(end);
        path
    } else {
        vec![
            start,
            (start.0 + 10.0, lane_y_mm),
            (end.0 - 10.0, lane_y_mm),
            end,
        ]
    };
    apply_series_channel_metadata(built, &channel.route, None)
}

pub fn create_parallel_geometry_from_spec(spec: &BlueprintTopologySpec) -> NetworkBlueprint {
    let inlet = layout_node("inlet", NodeKind::Inlet, (0.0, spec.box_dims_mm.1 * 0.5));
    let outlet = layout_node(
        "outlet",
        NodeKind::Outlet,
        (spec.box_dims_mm.0, spec.box_dims_mm.1 * 0.5),
    );
    let channel_count = spec.parallel_channels.len().max(1);
    let lane_spacing = spec.box_dims_mm.1 / (channel_count as f64 + 1.0);
    let channels = spec
        .parallel_channels
        .iter()
        .enumerate()
        .map(|(idx, channel)| {
            build_parallel_channel(channel, &inlet, &outlet, lane_spacing * (idx as f64 + 1.0))
        })
        .collect();

    let mut blueprint = NetworkBlueprint {
        name: spec.design_name.clone(),
        box_dims: spec.box_dims_mm,
        box_outline: make_outline(spec.box_dims_mm),
        nodes: vec![inlet, outlet],
        channels,
        render_hints: None,
        topology: Some(spec.clone()),
        lineage: None,
        metadata: None,
        geometry_authored: true,
    };
    blueprint.insert_metadata(GeometryAuthoringProvenance::create_geometry());
    materialize_blueprint_layout(&mut blueprint);
    blueprint
}
