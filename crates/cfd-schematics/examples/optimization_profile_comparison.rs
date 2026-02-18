use cfd_schematics::{
    config::{presets, ChannelTypeConfig, GeometryConfig, OptimizationProfile, SerpentineConfig},
    geometry::{generator::create_geometry, optimization::calculate_path_length, SplitType},
    visualizations::schematic::plot_geometry,
};
use std::fs;

fn main() {
    fs::create_dir_all("outputs/optimization").unwrap();

    let config = GeometryConfig::default();
    let splits = vec![SplitType::Bifurcation];
    let box_dims = (200.0, 100.0);

    println!("Serpentine Optimization Profile Comparison");
    println!("==========================================");
    println!();

    println!("1. Standard Configuration (No Optimization)");
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

    plot_geometry(
        &standard_system,
        "outputs/optimization/profile_standard.png",
    )
    .unwrap();
    println!("   Generation time: {:?}", standard_time);
    println!("   Total length:    {:.2} mm", standard_length);
    println!("   Generated: outputs/optimization/profile_standard.png");
    println!();

    println!("2. Fast Optimization Profile");
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

    plot_geometry(&fast_system, "outputs/optimization/profile_fast.png").unwrap();
    println!("   Generation time: {:?}", fast_time);
    println!("   Total length:    {:.2} mm", fast_length);
    println!("   Improvement:     {:.1}%", fast_improvement);
    println!(
        "   Speed ratio:     {:.1}x slower",
        fast_time.as_secs_f64() / standard_time.as_secs_f64()
    );
    println!("   Generated: outputs/optimization/profile_fast.png");
    println!();

    println!("3. Balanced Optimization Profile");
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

    plot_geometry(
        &balanced_system,
        "outputs/optimization/profile_balanced.png",
    )
    .unwrap();
    println!("   Generation time: {:?}", balanced_time);
    println!("   Total length:    {:.2} mm", balanced_length);
    println!("   Improvement:     {:.1}%", balanced_improvement);
    println!(
        "   Speed ratio:     {:.1}x slower",
        balanced_time.as_secs_f64() / standard_time.as_secs_f64()
    );
    println!("   Generated: outputs/optimization/profile_balanced.png");
    println!();

    println!("4. Thorough Optimization Profile");
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

    plot_geometry(
        &thorough_system,
        "outputs/optimization/profile_thorough.png",
    )
    .unwrap();
    println!("   Generation time: {:?}", thorough_time);
    println!("   Total length:    {:.2} mm", thorough_length);
    println!("   Improvement:     {:.1}%", thorough_improvement);
    println!(
        "   Speed ratio:     {:.1}x slower",
        thorough_time.as_secs_f64() / standard_time.as_secs_f64()
    );
    println!("   Generated: outputs/optimization/profile_thorough.png");
    println!();

    println!("5. Custom Configuration with Fast Profile");
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

    plot_geometry(&custom_system, "outputs/optimization/profile_custom.png").unwrap();
    println!("   Generation time: {:?}", custom_time);
    println!("   Total length:    {:.2} mm", custom_length);
    println!("   Improvement:     {:.1}%", custom_improvement);
    println!(
        "   Speed ratio:     {:.1}x slower",
        custom_time.as_secs_f64() / standard_time.as_secs_f64()
    );
    println!("   Generated: outputs/optimization/profile_custom.png");
    println!();

    println!("Summary");
    println!("=======");
    println!("Configuration        | Length (mm) | Improvement | Speed Ratio");
    println!("---------------------|-------------|-------------|------------");
    println!(
        "Standard (no opt)    | {:>9.2}   | {:>7.1}%   | {:>7.1}x",
        standard_length, 0.0, 1.0
    );
    println!(
        "Fast optimization    | {:>9.2}   | {:>7.1}%   | {:>7.1}x",
        fast_length,
        fast_improvement,
        fast_time.as_secs_f64() / standard_time.as_secs_f64()
    );
    println!(
        "Balanced optimization| {:>9.2}   | {:>7.1}%   | {:>7.1}x",
        balanced_length,
        balanced_improvement,
        balanced_time.as_secs_f64() / standard_time.as_secs_f64()
    );
    println!(
        "Thorough optimization| {:>9.2}   | {:>7.1}%   | {:>7.1}x",
        thorough_length,
        thorough_improvement,
        thorough_time.as_secs_f64() / standard_time.as_secs_f64()
    );
    println!(
        "Custom (fast profile)| {:>9.2}   | {:>7.1}%   | {:>7.1}x",
        custom_length,
        custom_improvement,
        custom_time.as_secs_f64() / standard_time.as_secs_f64()
    );
    println!();

    println!("Recommendations");
    println!("===============");
    println!("• Fast Profile:     Use for real-time applications where speed is critical");
    println!("• Balanced Profile: Good general-purpose optimization for most use cases");
    println!("• Thorough Profile: Use for final designs where maximum length is needed");
    println!("• Custom Config:    Fine-tune parameters for specific requirements");
    println!();
    println!("All profile comparison files generated in outputs/optimization/");
}

fn calculate_total_length(system: &cfd_schematics::geometry::ChannelSystem) -> f64 {
    system
        .channels
        .iter()
        .map(|channel| match &channel.channel_type {
            cfd_schematics::geometry::ChannelType::Serpentine { path } => {
                calculate_path_length(path)
            }
            _ => 0.0,
        })
        .sum()
}
