use cfd_schematics::{
    config::{presets, ChannelTypeConfig, GeometryConfig, OptimizationProfile, SerpentineConfig},
    geometry::{generator::create_geometry, optimization::calculate_path_length, SplitType},
};
#[path = "../shared/mod.rs"]
mod shared;

fn main() {

    let config = GeometryConfig::default();
    let splits = vec![SplitType::Bifurcation];
    let box_dims = (200.0, 100.0);

    tracing::info!("Serpentine Optimization Profile Comparison");
    tracing::info!("==========================================");

    tracing::info!("1. Standard Configuration (No Optimization)");
    let start_time = std::time::Instant::now();
    let standard_config = SerpentineConfig::default();
    let standard_system = create_geometry(
        box_dims,
        &splits,
        &config,
        &ChannelTypeConfig::AllSerpentine(standard_config),
    );
    let standard_time = start_time.elapsed();
    let standard_length = calculate_total_length(&standard_system);

    shared::save_example_output_with_name(&standard_system, "optimization_profile_comparison", "profile_standard");

    tracing::info!("   Generation time: {:?}", standard_time);
    tracing::info!("   Total length:    {:.2} mm", standard_length);

    tracing::info!("2. Fast Optimization Profile");
    let start_time = std::time::Instant::now();
    let fast_config = presets::fast_optimized_serpentine();
    let fast_system = create_geometry(
        box_dims,
        &splits,
        &config,
        &ChannelTypeConfig::AllSerpentine(fast_config),
    );
    let fast_time = start_time.elapsed();
    let fast_length = calculate_total_length(&fast_system);
    let fast_improvement = ((fast_length - standard_length) / standard_length) * 100.0;

    shared::save_example_output_with_name(&fast_system, "optimization_profile_comparison", "profile_fast");

    tracing::info!("   Generation time: {:?}", fast_time);
    tracing::info!("   Total length:    {:.2} mm", fast_length);
    tracing::info!("   Improvement:     {:.1}%", fast_improvement);
    tracing::info!(
        "   Speed ratio:     {:.1}x slower",
        fast_time.as_secs_f64() / standard_time.as_secs_f64()
    );

    tracing::info!("3. Balanced Optimization Profile");
    let start_time = std::time::Instant::now();
    let balanced_config = presets::optimized_serpentine();
    let balanced_system = create_geometry(
        box_dims,
        &splits,
        &config,
        &ChannelTypeConfig::AllSerpentine(balanced_config),
    );
    let balanced_time = start_time.elapsed();
    let balanced_length = calculate_total_length(&balanced_system);
    let balanced_improvement = ((balanced_length - standard_length) / standard_length) * 100.0;

    shared::save_example_output_with_name(&balanced_system, "optimization_profile_comparison", "profile_balanced");

    tracing::info!("   Generation time: {:?}", balanced_time);
    tracing::info!("   Total length:    {:.2} mm", balanced_length);
    tracing::info!("   Improvement:     {:.1}%", balanced_improvement);
    tracing::info!(
        "   Speed ratio:     {:.1}x slower",
        balanced_time.as_secs_f64() / standard_time.as_secs_f64()
    );

    tracing::info!("4. Thorough Optimization Profile");
    let start_time = std::time::Instant::now();
    let thorough_config = presets::thorough_optimized_serpentine();
    let thorough_system = create_geometry(
        box_dims,
        &splits,
        &config,
        &ChannelTypeConfig::AllSerpentine(thorough_config),
    );
    let thorough_time = start_time.elapsed();
    let thorough_length = calculate_total_length(&thorough_system);
    let thorough_improvement = ((thorough_length - standard_length) / standard_length) * 100.0;

    shared::save_example_output_with_name(&thorough_system, "optimization_profile_comparison", "profile_thorough");

    tracing::info!("   Generation time: {:?}", thorough_time);
    tracing::info!("   Total length:    {:.2} mm", thorough_length);
    tracing::info!("   Improvement:     {:.1}%", thorough_improvement);
    tracing::info!(
        "   Speed ratio:     {:.1}x slower",
        thorough_time.as_secs_f64() / standard_time.as_secs_f64()
    );

    tracing::info!("5. Custom Configuration with Fast Profile");
    let start_time = std::time::Instant::now();
    let custom_config = SerpentineConfig {
        fill_factor: 0.85,
        wavelength_factor: 2.5,
        gaussian_width_factor: 8.0,
        wave_density_factor: 3.5,
        wave_phase_direction: 0.0,
        optimization_enabled: true,
        target_fill_ratio: 0.92,
        optimization_profile: OptimizationProfile::Fast,
        ..SerpentineConfig::default()
    };
    let custom_system = create_geometry(
        box_dims,
        &splits,
        &config,
        &ChannelTypeConfig::AllSerpentine(custom_config),
    );
    let custom_time = start_time.elapsed();
    let custom_length = calculate_total_length(&custom_system);
    let custom_improvement = ((custom_length - standard_length) / standard_length) * 100.0;

    shared::save_example_output_with_name(&custom_system, "optimization_profile_comparison", "profile_custom");

    tracing::info!("   Generation time: {:?}", custom_time);
    tracing::info!("   Total length:    {:.2} mm", custom_length);
    tracing::info!("   Improvement:     {:.1}%", custom_improvement);
    tracing::info!(
        "   Speed ratio:     {:.1}x slower",
        custom_time.as_secs_f64() / standard_time.as_secs_f64()
    );

    tracing::info!("Summary");
    tracing::info!("=======");
    tracing::info!("Configuration        | Length (mm) | Improvement | Speed Ratio");
    tracing::info!("---------------------|-------------|-------------|------------");
    tracing::info!(
        "Standard (no opt)    | {:>9.2}   | {:>7.1}%   | {:>7.1}x",
        standard_length,
        0.0,
        1.0
    );
    tracing::info!(
        "Fast optimization    | {:>9.2}   | {:>7.1}%   | {:>7.1}x",
        fast_length,
        fast_improvement,
        fast_time.as_secs_f64() / standard_time.as_secs_f64()
    );
    tracing::info!(
        "Balanced optimization| {:>9.2}   | {:>7.1}%   | {:>7.1}x",
        balanced_length,
        balanced_improvement,
        balanced_time.as_secs_f64() / standard_time.as_secs_f64()
    );
    tracing::info!(
        "Thorough optimization| {:>9.2}   | {:>7.1}%   | {:>7.1}x",
        thorough_length,
        thorough_improvement,
        thorough_time.as_secs_f64() / standard_time.as_secs_f64()
    );
    tracing::info!(
        "Custom (fast profile)| {:>9.2}   | {:>7.1}%   | {:>7.1}x",
        custom_length,
        custom_improvement,
        custom_time.as_secs_f64() / standard_time.as_secs_f64()
    );

    tracing::info!("Recommendations");
    tracing::info!("===============");
    tracing::info!("• Fast Profile:     Use for real-time applications where speed is critical");
    tracing::info!("• Balanced Profile: Good general-purpose optimization for most use cases");
    tracing::info!("• Thorough Profile: Use for final designs where maximum length is needed");
    tracing::info!("• Custom Config:    Fine-tune parameters for specific requirements");
    tracing::info!("All profile comparison files generated in output/examples/optimization_profile_comparison/");
}

fn calculate_total_length(system: &cfd_schematics::NetworkBlueprint) -> f64 {
    system
        .channels
        .iter()
        .map(|channel| match &channel.channel_shape {
            cfd_schematics::domain::model::ChannelShape::Serpentine { .. } => {
                calculate_path_length(&channel.path)
            }
            _ => 0.0,
        })
        .sum()
}
