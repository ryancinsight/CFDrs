//! Compatibility with previously serialized simulation metadata.

use super::SimulationMetadata;
use serde_json::json;
use tracing_subscriber::fmt::format::FmtSpan;

#[test]
fn timestamp_fields_keep_the_utc_json_contract() {
    let _trace = tracing::subscriber::set_default(
        tracing_subscriber::fmt()
            .with_span_events(FmtSpan::FULL)
            .with_test_writer()
            .finish(),
    );
    let prior = tracing::info_span!("arrange").in_scope(|| {
        json!({
            "name": "Flow", "description": "", "version": "1.0.0", "author": "",
            "tags": [], "created_at": "2020-02-29T01:02:03.250Z",
            "modified_at": "2020-03-01T01:02:03-05:00"
        })
    });
    let metadata: SimulationMetadata = tracing::info_span!("act")
        .in_scope(|| serde_json::from_value(prior).expect("parse existing simulation metadata"));
    let _assert = tracing::info_span!("assert").entered();
    assert!(metadata.modified_at > metadata.created_at);
    let current = serde_json::to_value(metadata).expect("serialize simulation metadata");
    assert_eq!(current["created_at"], "2020-02-29T01:02:03.250Z");
    assert_eq!(current["modified_at"], "2020-03-01T06:02:03Z");
}
