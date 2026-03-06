use crate::domain::model::NetworkBlueprint;
use crate::geometry::metadata::{
    ChannelPathMetadata, ChannelVisualRole, JunctionFamily, JunctionGeometryMetadata,
    NodeLayoutMetadata,
};

const PLATE_W_MM: f64 = 127.76;
const PLATE_H_MM: f64 = 85.47;
const Y_MID_MM: f64 = PLATE_H_MM * 0.5;

pub(crate) fn attach_asymmetric_trifurcation_layout(bp: &mut NetworkBlueprint) {
    set_node(bp, "inlet", 8.0, Y_MID_MM);
    set_node(bp, "split_jn", 26.0, Y_MID_MM);
    set_node(bp, "throat_in", 62.0, Y_MID_MM);
    set_node(bp, "throat_out", 68.0, Y_MID_MM);
    set_node(bp, "periph_merge", 98.0, Y_MID_MM);
    set_node(bp, "outlet_merge", 112.0, Y_MID_MM);
    set_node(bp, "outlet", 124.0, Y_MID_MM);
    set_junction(bp, "split_jn", JunctionFamily::Trifurcation, &[-25.0, 0.0, 30.0], &[]);
    set_junction(bp, "periph_merge", JunctionFamily::Merge, &[], &[30.0, 25.0]);
    set_junction(bp, "outlet_merge", JunctionFamily::Merge, &[], &[18.0, 18.0]);
    set_polyline(
        bp,
        "left_arm",
        vec![
            node_point(bp, "split_jn"),
            (70.0, Y_MID_MM - 14.0),
            node_point(bp, "periph_merge"),
        ],
        ChannelVisualRole::PeripheralBypass,
    );
    set_polyline(
        bp,
        "right_arm",
        vec![
            node_point(bp, "split_jn"),
            (74.0, Y_MID_MM + 16.0),
            node_point(bp, "periph_merge"),
        ],
        ChannelVisualRole::PeripheralBypass,
    );
}

fn set_node(bp: &mut NetworkBlueprint, node_id: &str, x_mm: f64, y_mm: f64) {
    if let Some(node) = bp.nodes.iter_mut().find(|node| node.id.as_str() == node_id) {
        node.layout = Some(NodeLayoutMetadata { x_mm, y_mm });
        node.metadata
            .get_or_insert_with(crate::geometry::metadata::MetadataContainer::new)
            .insert(NodeLayoutMetadata { x_mm, y_mm });
    }
}

fn set_junction(
    bp: &mut NetworkBlueprint,
    node_id: &str,
    junction_family: JunctionFamily,
    branch_angles_deg: &[f64],
    merge_angles_deg: &[f64],
) {
    if let Some(node) = bp.nodes.iter_mut().find(|node| node.id.as_str() == node_id) {
        let geometry = JunctionGeometryMetadata {
            junction_family,
            branch_angles_deg: branch_angles_deg.to_vec(),
            merge_angles_deg: merge_angles_deg.to_vec(),
        };
        node.junction_geometry = Some(geometry.clone());
        node.metadata
            .get_or_insert_with(crate::geometry::metadata::MetadataContainer::new)
            .insert(geometry);
    }
}

fn set_polyline(
    bp: &mut NetworkBlueprint,
    channel_id: &str,
    polyline_mm: Vec<(f64, f64)>,
    visual_role: ChannelVisualRole,
) {
    if let Some(channel) = bp
        .channels
        .iter_mut()
        .find(|channel| channel.id.as_str() == channel_id)
    {
        let path = ChannelPathMetadata {
            polyline_mm: polyline_mm.clone(),
            visual_role,
        };
        channel.path = Some(path.clone());
        channel
            .metadata
            .get_or_insert_with(crate::geometry::metadata::MetadataContainer::new)
            .insert(path);
    }
}

fn node_point(bp: &NetworkBlueprint, node_id: &str) -> (f64, f64) {
    bp.nodes
        .iter()
        .find(|node| node.id.as_str() == node_id)
        .and_then(|node| node.layout.map(|layout| (layout.x_mm, layout.y_mm)))
        .unwrap_or((PLATE_W_MM * 0.5, Y_MID_MM))
}
