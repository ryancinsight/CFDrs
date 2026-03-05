//! Lightweight strict template renderer for Milestone 12 narrative markdown.
//!
//! Placeholders use `{{TOKEN}}` form. Rendering fails if any placeholder
//! remains unresolved to prevent silent drift in milestone reports.

use std::collections::{BTreeMap, BTreeSet};

/// Render a template using strict placeholder replacement.
///
/// # Errors
/// Returns an error if unresolved placeholders remain in the rendered output.
pub fn render_template_strict(
    template: &str,
    values: &BTreeMap<String, String>,
) -> Result<String, Box<dyn std::error::Error>> {
    let mut rendered = template.to_string();
    for (key, value) in values {
        let token = format!("{{{{{key}}}}}");
        rendered = rendered.replace(&token, value);
    }

    let unresolved = extract_placeholders(&rendered);
    if unresolved.is_empty() {
        Ok(rendered)
    } else {
        Err(format!(
            "unresolved template placeholders: {}",
            unresolved.into_iter().collect::<Vec<_>>().join(", ")
        )
        .into())
    }
}

/// Extract all `{{TOKEN}}` placeholders from a template string.
#[must_use]
pub fn extract_placeholders(template: &str) -> BTreeSet<String> {
    let mut out = BTreeSet::new();
    let mut idx = 0usize;
    while let Some(open_rel) = template[idx..].find("{{") {
        let open = idx + open_rel;
        let search_start = open + 2;
        if let Some(close_rel) = template[search_start..].find("}}") {
            let close = search_start + close_rel;
            let token = template[search_start..close].trim();
            if !token.is_empty() {
                out.insert(token.to_string());
            }
            idx = close + 2;
        } else {
            break;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strict_renderer_replaces_all_tokens() {
        let tpl = "A {{ONE}} B {{TWO}}";
        let mut map = BTreeMap::new();
        map.insert("ONE".to_string(), "1".to_string());
        map.insert("TWO".to_string(), "2".to_string());

        let rendered = render_template_strict(tpl, &map).expect("template should render");
        assert_eq!(rendered, "A 1 B 2");
    }

    #[test]
    fn strict_renderer_fails_on_unresolved() {
        let tpl = "A {{ONE}} B {{TWO}}";
        let mut map = BTreeMap::new();
        map.insert("ONE".to_string(), "1".to_string());
        let err = render_template_strict(tpl, &map).expect_err("must fail");
        assert!(err.to_string().contains("TWO"));
    }
}
