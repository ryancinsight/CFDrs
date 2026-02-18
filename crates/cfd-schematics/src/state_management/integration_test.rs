//! Integration tests for the state management system
//!
//! This module provides basic integration tests to verify that the
//! parameter registry and managers work correctly together.

#[cfg(test)]
mod tests {
    use super::super::{adaptive::ChannelGenerationContext, ParameterManager, ParameterRegistry};
    use crate::config::GeometryConfig;

    #[test]
    fn test_parameter_registry_creation() {
        let registry = ParameterRegistry::with_defaults();
        assert!(
            registry.is_ok(),
            "Failed to create parameter registry: {:?}",
            registry.err()
        );

        let registry = registry.unwrap();
        assert_eq!(
            registry.domain_names(),
            vec!["serpentine", "arc", "geometry", "collision", "symmetry"]
        );
        assert!(registry.is_validation_enabled());
        assert!(!registry.is_locked());
    }

    #[test]
    fn test_serpentine_parameter_manager() {
        let registry = ParameterRegistry::with_defaults().unwrap();
        let serpentine_manager = registry.serpentine();

        // Test parameter names
        let param_names = serpentine_manager.parameter_names();
        assert!(param_names.contains(&"amplitude".to_string()));
        assert!(param_names.contains(&"wavelength_factor".to_string()));
        assert!(param_names.contains(&"wave_density_factor".to_string()));

        // Test parameter existence
        assert!(serpentine_manager.has_parameter("amplitude"));
        assert!(serpentine_manager.has_parameter("wavelength_factor"));
        assert!(!serpentine_manager.has_parameter("nonexistent_param"));

        // Test validation
        assert!(serpentine_manager.validate_all().is_ok());
    }

    #[test]
    fn test_arc_parameter_manager() {
        let registry = ParameterRegistry::with_defaults().unwrap();
        let arc_manager = registry.arc();

        // Test parameter names
        let param_names = arc_manager.parameter_names();
        assert!(param_names.contains(&"curvature_factor".to_string()));
        assert!(param_names.contains(&"smoothness".to_string()));

        // Test parameter existence
        assert!(arc_manager.has_parameter("curvature_factor"));
        assert!(arc_manager.has_parameter("smoothness"));
        assert!(!arc_manager.has_parameter("nonexistent_param"));

        // Test validation
        assert!(arc_manager.validate_all().is_ok());
    }

    #[test]
    fn test_geometry_parameter_manager() {
        let registry = ParameterRegistry::with_defaults().unwrap();
        let geometry_manager = registry.geometry();

        // Test parameter names
        let param_names = geometry_manager.parameter_names();
        assert!(param_names.contains(&"wall_clearance".to_string()));
        assert!(param_names.contains(&"channel_width".to_string()));
        assert!(param_names.contains(&"channel_height".to_string()));

        // Test parameter existence
        assert!(geometry_manager.has_parameter("wall_clearance"));
        assert!(geometry_manager.has_parameter("channel_width"));
        assert!(!geometry_manager.has_parameter("nonexistent_param"));

        // Test validation
        assert!(geometry_manager.validate_all().is_ok());
    }

    #[test]
    fn test_adaptive_parameter_behavior() {
        let registry = ParameterRegistry::with_defaults().unwrap();
        let serpentine_manager = registry.serpentine();

        // Create a test context
        let context = ChannelGenerationContext::new(
            GeometryConfig::default(),
            (100.0, 50.0),
            4,
            Some(&[10.0, 20.0, 30.0, 40.0]),
        )
        .with_endpoints((0.0, 25.0), (100.0, 25.0));

        // Test adaptive amplitude
        let base_amplitude = serpentine_manager.get_amplitude(None);
        let adaptive_amplitude = serpentine_manager.get_amplitude(Some(&context));

        // The adaptive amplitude should be different (likely smaller due to constraints)
        // but both should be positive
        assert!(base_amplitude > 0.0);
        assert!(adaptive_amplitude > 0.0);

        // Test adaptive wavelength factor
        let base_wavelength = serpentine_manager.get_wavelength_factor(None);
        let adaptive_wavelength = serpentine_manager.get_wavelength_factor(Some(&context));

        assert!(base_wavelength > 0.0);
        assert!(adaptive_wavelength > 0.0);
    }

    #[test]
    fn test_parameter_metadata() {
        let registry = ParameterRegistry::with_defaults().unwrap();
        let serpentine_manager = registry.serpentine();

        // Test amplitude metadata
        let amplitude_metadata = serpentine_manager.get_metadata("amplitude");
        assert!(amplitude_metadata.is_ok());

        let metadata = amplitude_metadata.unwrap();
        assert_eq!(metadata.name, "amplitude");
        assert_eq!(metadata.category, "wave_parameters");
        assert!(metadata.is_mutable);
        assert!(metadata.affects_others);
        assert_eq!(metadata.units, Some("mm".to_string()));
    }

    #[test]
    fn test_registry_validation() {
        let registry = ParameterRegistry::with_defaults().unwrap();

        // Test global validation
        assert!(registry.validate_all().is_ok());

        // Test that validation is enabled by default
        assert!(registry.is_validation_enabled());
    }

    #[test]
    fn test_registry_locking() {
        let mut registry = ParameterRegistry::with_defaults().unwrap();

        // Initially unlocked
        assert!(!registry.is_locked());

        // Lock the registry
        registry.lock();
        assert!(registry.is_locked());

        // Trying to get mutable access should fail
        assert!(registry.serpentine_mut().is_err());

        // Unlock the registry
        registry.unlock();
        assert!(!registry.is_locked());

        // Now mutable access should work
        assert!(registry.serpentine_mut().is_ok());
    }

    #[test]
    fn test_parameter_names_across_domains() {
        let registry = ParameterRegistry::with_defaults().unwrap();
        let all_params = registry.all_parameter_names();

        // Check that we have parameters for each domain
        assert!(all_params.contains_key("serpentine"));
        assert!(all_params.contains_key("arc"));
        assert!(all_params.contains_key("geometry"));

        // Check that serpentine has expected parameters
        let serpentine_params = &all_params["serpentine"];
        assert!(serpentine_params.contains(&"amplitude".to_string()));
        assert!(serpentine_params.contains(&"wavelength_factor".to_string()));

        // Check that arc has expected parameters
        let arc_params = &all_params["arc"];
        assert!(arc_params.contains(&"curvature_factor".to_string()));
        assert!(arc_params.contains(&"smoothness".to_string()));

        // Check that geometry has expected parameters
        let geometry_params = &all_params["geometry"];
        assert!(geometry_params.contains(&"wall_clearance".to_string()));
        assert!(geometry_params.contains(&"channel_width".to_string()));
    }

    #[test]
    fn test_wave_parameters_helper() {
        let registry = ParameterRegistry::with_defaults().unwrap();
        let serpentine_manager = registry.serpentine();

        // Test the wave parameters helper method
        let wave_params = serpentine_manager.get_wave_parameters(None);

        // Check that all expected parameters are present
        assert!(wave_params.contains_key("amplitude"));
        assert!(wave_params.contains_key("wavelength_factor"));
        assert!(wave_params.contains_key("frequency_multiplier"));
        assert!(wave_params.contains_key("phase_offset"));
        assert!(wave_params.contains_key("gaussian_width_factor"));
        assert!(wave_params.contains_key("wave_density_factor"));
        assert!(wave_params.contains_key("fill_factor"));
        assert!(wave_params.contains_key("target_fill_ratio"));

        // Check that all values are reasonable
        for (name, value) in wave_params {
            assert!(
                value.is_finite(),
                "Parameter {} has invalid value: {}",
                name,
                value
            );
            if name != "phase_offset" {
                // phase_offset can be 0
                assert!(
                    value > 0.0,
                    "Parameter {} should be positive: {}",
                    name,
                    value
                );
            }
        }
    }
}
