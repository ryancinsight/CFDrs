//! Plugin storage management.

use crate::error::{Error, PluginErrorKind, Result};
use crate::plugin::traits::Plugin;
use std::collections::HashMap;
use std::sync::Arc;

/// Standard data storage for plugins (no internal locking)
pub(crate) struct PluginStorage {
    plugins: HashMap<String, Arc<dyn Plugin>>,
}

impl PluginStorage {
    pub fn new() -> Self {
        Self {
            plugins: HashMap::new(),
        }
    }

    pub fn insert(&mut self, plugin: Arc<dyn Plugin>) -> Result<()> {
        let name = plugin.name().to_string();
        if self.plugins.contains_key(&name) {
            return Err(Error::Plugin(PluginErrorKind::AlreadyRegistered { name }));
        }
        self.plugins.insert(name, plugin);
        Ok(())
    }

    pub fn get(&self, name: &str) -> Option<&Arc<dyn Plugin>> {
        self.plugins.get(name)
    }

    pub fn remove(&mut self, name: &str) -> Option<Arc<dyn Plugin>> {
        self.plugins.remove(name)
    }

    pub fn list(&self) -> Vec<String> {
        self.plugins.keys().cloned().collect()
    }

    pub fn contains(&self, name: &str) -> bool {
        self.plugins.contains_key(name)
    }

    pub fn iter(&self) -> impl Iterator<Item = (&String, &Arc<dyn Plugin>)> {
        self.plugins.iter()
    }

    pub fn clear(&mut self) {
        self.plugins.clear();
    }

    pub fn len(&self) -> usize {
        self.plugins.len()
    }

    pub fn is_empty(&self) -> bool {
        self.plugins.is_empty()
    }
}
