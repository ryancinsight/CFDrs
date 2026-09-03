use crate::domain::model::{ChannelSpec, NetworkBlueprint};
use crate::domain::therapy_metadata::TherapyZone;
use crate::geometry::metadata::ChannelVisualRole;
use crate::geometry::Point2D;
use std::borrow::Cow;
use std::collections::BTreeMap;

pub(super) fn treatment_lane_paths<'blueprint>(
    blueprint: &'blueprint NetworkBlueprint,
) -> Vec<(f64, Cow<'blueprint, [Point2D]>)> {
    type LaneSegment<'a> = (f64, Cow<'a, [Point2D]>);

    let mut grouped: BTreeMap<i32, Vec<LaneSegment<'blueprint>>> = BTreeMap::new();
    for channel in &blueprint.channels {
        if channel.therapy_zone != Some(TherapyZone::CancerTarget)
            || channel.visual_role == Some(ChannelVisualRole::Trunk)
            || channel.path.len() < 2
        {
            continue;
        }
        let centroid_y = channel_centroid_y(channel);
        let bucket = (centroid_y * 2.0).round() as i32;
        let (x_min, _) = x_span(&channel.path);
        let path = match (channel.path.first().copied(), channel.path.last().copied()) {
            (Some(first), Some(last)) if first.0 > last.0 => {
                let mut reversed = channel.path.clone();
                reversed.reverse();
                Cow::Owned(reversed)
            }
            _ => Cow::Borrowed(channel.path.as_slice()),
        };
        grouped.entry(bucket).or_default().push((x_min, path));
    }

    grouped
        .into_iter()
        .map(|(bucket, mut segments)| {
            segments.sort_by(|left, right| left.0.total_cmp(&right.0));
            if segments.len() == 1 {
                let Some((_, segment)) = segments.pop() else {
                    return (f64::from(bucket) / 2.0, Cow::Owned(Vec::new()));
                };
                return (f64::from(bucket) / 2.0, segment);
            }

            let lane_capacity = segments.iter().fold(0usize, |capacity, (_, segment)| {
                capacity.saturating_add(segment.len())
            });
            let mut lane = Vec::with_capacity(lane_capacity);
            for (_, segment) in segments {
                if lane
                    .last()
                    .zip(segment.as_ref().first())
                    .is_some_and(|(last, first)| *last == *first)
                {
                    lane.extend(segment.as_ref().iter().skip(1).copied());
                } else {
                    lane.extend(segment.as_ref().iter().copied());
                }
            }
            (f64::from(bucket) / 2.0, Cow::Owned(lane))
        })
        .filter(|(_, lane)| lane.len() >= 2)
        .collect()
}

pub(super) fn channel_centroid_y(channel: &ChannelSpec) -> f64 {
    channel.path.iter().map(|(_, y)| *y).sum::<f64>() / channel.path.len() as f64
}

pub(super) fn x_span(path: &[Point2D]) -> (f64, f64) {
    path.iter().fold(
        (f64::INFINITY, f64::NEG_INFINITY),
        |(x_min, x_max), (x, _)| (x_min.min(*x), x_max.max(*x)),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn treatment_channel(path: Vec<Point2D>) -> ChannelSpec {
        let mut channel =
            ChannelSpec::new_pipe_rect("treatment", "inlet", "outlet", 1.0, 1.0, 1.0, 1.0, 0.0)
                .with_path(path);
        channel.therapy_zone = Some(TherapyZone::CancerTarget);
        channel
    }

    #[test]
    fn single_forward_treatment_lane_borrows_channel_path() {
        let path = vec![(0.0, 10.0), (100.0, 10.0)];
        let mut blueprint = NetworkBlueprint::new_with_explicit_positions("lane");
        blueprint.channels.push(treatment_channel(path.clone()));

        let lanes = treatment_lane_paths(&blueprint);
        assert_eq!(lanes.len(), 1);
        let Cow::Borrowed(lane_path) = &lanes[0].1 else {
            panic!("single forward lane must remain borrowed");
        };
        assert_eq!(*lane_path, path.as_slice());
    }

    #[test]
    fn reversed_treatment_lane_is_owned_in_forward_order() {
        let mut blueprint = NetworkBlueprint::new_with_explicit_positions("lane");
        blueprint
            .channels
            .push(treatment_channel(vec![(100.0, 10.0), (0.0, 10.0)]));

        let lanes = treatment_lane_paths(&blueprint);
        assert_eq!(lanes.len(), 1);
        assert!(matches!(&lanes[0].1, Cow::Owned(_)));
        assert_eq!(lanes[0].1.as_ref(), &[(0.0, 10.0), (100.0, 10.0)]);
    }

    #[test]
    fn merged_treatment_lane_preserves_segment_order() {
        let mut blueprint = NetworkBlueprint::new_with_explicit_positions("lane");
        blueprint
            .channels
            .push(treatment_channel(vec![(0.0, 10.0), (50.0, 10.0)]));
        blueprint
            .channels
            .push(treatment_channel(vec![(50.0, 10.0), (100.0, 10.0)]));

        let lanes = treatment_lane_paths(&blueprint);
        assert_eq!(lanes.len(), 1);
        assert!(matches!(&lanes[0].1, Cow::Owned(_)));
        assert_eq!(
            lanes[0].1.as_ref(),
            &[(0.0, 10.0), (50.0, 10.0), (100.0, 10.0)]
        );
    }
}
