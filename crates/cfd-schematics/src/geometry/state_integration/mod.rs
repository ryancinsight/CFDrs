//! Integration layer between existing strategies and new state management
//!
//! This module provides integration utilities to gradually migrate from
//! hardcoded parameters to the centralized state management system while
//! maintaining backward compatibility.

mod arc;
mod serpentine;

pub use arc::*;
pub use serpentine::*;
