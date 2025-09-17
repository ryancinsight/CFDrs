//! Dependency resolution for plugins.

use crate::error::{Error, PluginErrorKind, Result};
use std::collections::{HashMap, HashSet};

/// Standard dependency resolver (no internal locking)
pub(crate) struct DependencyResolver {
    dependencies: HashMap<String, Vec<String>>,
    load_order: Vec<String>,
}

impl DependencyResolver {
    pub fn new() -> Self {
        Self {
            dependencies: HashMap::new(),
            load_order: Vec::new(),
        }
    }

    pub fn add_plugin(&mut self, name: String, deps: Vec<String>) -> Result<()> {
        self.dependencies.insert(name, deps);
        self.update_load_order()
    }

    pub fn remove_plugin(&mut self, name: &str) {
        self.dependencies.remove(name);
        let _ = self.update_load_order(); // Ignore errors during removal
    }

    pub fn validate_dependencies(
        &self,
        deps: &[String],
        available_plugins: &[String],
    ) -> Result<()> {
        for dep in deps {
            if !available_plugins.contains(dep) {
                return Err(Error::Plugin(PluginErrorKind::DependencyMissing {
                    plugin: String::new(),
                    dependency: dep.to_string(),
                }));
            }
        }
        Ok(())
    }

    pub fn get_load_order(&self) -> &[String] {
        &self.load_order
    }

    fn update_load_order(&mut self) -> Result<()> {
        // Detect circular dependencies first
        self.detect_circular_dependencies()?;

        // Perform topological sort
        let sorted = self.topological_sort()?;
        self.load_order = sorted;
        Ok(())
    }

    fn detect_circular_dependencies(&self) -> Result<()> {
        let mut visited = HashSet::new();
        let mut rec_stack = HashSet::new();

        for plugin in self.dependencies.keys() {
            if !visited.contains(plugin) && self.has_cycle(plugin, &mut visited, &mut rec_stack)? {
                return Err(Error::Plugin(PluginErrorKind::CircularDependency {
                    chain: self.find_cycle_chain(plugin),
                }));
            }
        }
        Ok(())
    }

    fn has_cycle(
        &self,
        node: &str,
        visited: &mut HashSet<String>,
        rec_stack: &mut HashSet<String>,
    ) -> Result<bool> {
        visited.insert(node.to_string());
        rec_stack.insert(node.to_string());

        if let Some(deps) = self.dependencies.get(node) {
            for dep in deps {
                if !visited.contains(dep) {
                    if self.has_cycle(dep, visited, rec_stack)? {
                        return Ok(true);
                    }
                } else if rec_stack.contains(dep) {
                    return Ok(true);
                }
            }
        }

        rec_stack.remove(node);
        Ok(false)
    }

    fn find_cycle_chain(&self, start: &str) -> Vec<String> {
        let mut chain = vec![start.to_string()];
        let mut current = start;

        while let Some(deps) = self.dependencies.get(current) {
            for dep in deps {
                if dep == start {
                    chain.push(dep.to_string());
                    return chain;
                }
                if !chain.contains(dep) {
                    chain.push(dep.to_string());
                    current = dep;
                    break;
                }
            }
        }
        chain
    }

    fn topological_sort(&self) -> Result<Vec<String>> {
        let mut in_degree: HashMap<String, usize> = HashMap::new();
        let mut adj_list: HashMap<String, Vec<String>> = HashMap::new();

        // Initialize in-degree and adjacency list
        for (plugin, deps) in &self.dependencies {
            in_degree.entry(plugin.to_string()).or_insert(0);
            for dep in deps {
                in_degree.entry(dep.to_string()).or_insert(0);
                adj_list
                    .entry(dep.to_string())
                    .or_default()
                    .push(plugin.to_string());
                *in_degree.get_mut(plugin).unwrap() += 1;
            }
        }

        // Find all nodes with in-degree 0
        let mut queue: Vec<String> = in_degree
            .iter()
            .filter_map(|(k, &v)| if v == 0 { Some(k.to_string()) } else { None })
            .collect();

        let mut result = Vec::new();

        while let Some(node) = queue.pop() {
            if let Some(neighbors) = adj_list.get(&node) {
                for neighbor in neighbors {
                    if let Some(degree) = in_degree.get_mut(neighbor) {
                        *degree -= 1;
                        if *degree == 0 {
                            queue.push(neighbor.to_string());
                        }
                    }
                }
            }
            result.push(node); // Move node instead of cloning
        }

        if result.len() != in_degree.len() {
            return Err(Error::Plugin(PluginErrorKind::CircularDependency {
                chain: vec!["Unable to determine exact cycle".to_string()],
            }));
        }

        Ok(result)
    }
}
