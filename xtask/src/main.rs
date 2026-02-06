//! xtask - Build automation for CFD-rs Python bindings and validation
//!
//! This tool automates the complete workflow:
//! 1. Build pycfdrs wheels with maturin
//! 2. Setup Python virtual environment
//! 3. Install pycfdrs and validation dependencies
//! 4. Install FEniCS or alternative CFD packages
//! 5. Run validation suite
//!
//! Usage:
//!   cargo xtask build-wheel       # Build pycfdrs wheel
//!   cargo xtask setup-venv        # Create virtual environment
//!   cargo xtask install-deps      # Install Python dependencies
//!   cargo xtask install-fenics    # Install FEniCS for validation
//!   cargo xtask validate          # Run validation suite
//!   cargo xtask all               # Complete workflow

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use std::env;
use std::path::{Path, PathBuf};
use xshell::{cmd, Shell};

#[derive(Parser)]
#[command(name = "xtask")]
#[command(about = "Build automation for CFD-rs validation")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Build pycfdrs wheel with maturin
    BuildWheel {
        /// Build in release mode
        #[arg(long, default_value_t = true)]
        release: bool,
    },

    /// Setup Python virtual environment
    SetupVenv {
        /// Python executable to use
        #[arg(long, default_value = "python")]
        python: String,

        /// Virtual environment path
        #[arg(long, default_value = ".venv")]
        path: PathBuf,
    },

    /// Install Python dependencies
    InstallDeps {
        /// Virtual environment path
        #[arg(long, default_value = ".venv")]
        venv: PathBuf,

        /// Include FEniCS
        #[arg(long)]
        with_fenics: bool,
    },

    /// Install FEniCS via conda (recommended)
    InstallFenics {
        /// Use conda instead of pip
        #[arg(long, default_value_t = true)]
        use_conda: bool,
    },

    /// Run validation suite
    Validate {
        /// Virtual environment path
        #[arg(long, default_value = ".venv")]
        venv: PathBuf,

        /// Generate plots
        #[arg(long)]
        plot: bool,

        /// Quick mode (coarse grids)
        #[arg(long)]
        quick: bool,

        /// Test category (analytical, fenics, literature, all)
        #[arg(long, default_value = "all")]
        category: String,
    },

    /// Complete workflow: build, setup, install, validate
    All {
        /// Use conda for FEniCS
        #[arg(long)]
        with_fenics: bool,

        /// Generate plots
        #[arg(long)]
        plot: bool,
    },

    /// Clean build artifacts and virtual environment
    Clean,

    /// Check if compilation errors exist
    Check,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let sh = Shell::new()?;

    match cli.command {
        Commands::BuildWheel { release } => build_wheel(&sh, release),
        Commands::SetupVenv { python, path } => setup_venv(&sh, &python, &path),
        Commands::InstallDeps { venv, with_fenics } => install_deps(&sh, &venv, with_fenics),
        Commands::InstallFenics { use_conda } => install_fenics(&sh, use_conda),
        Commands::Validate {
            venv,
            plot,
            quick,
            category,
        } => validate(&sh, &venv, plot, quick, &category),
        Commands::All { with_fenics, plot } => run_all(&sh, with_fenics, plot),
        Commands::Clean => clean(&sh),
        Commands::Check => check_compilation(&sh),
    }
}

/// Get project root directory
fn project_root() -> PathBuf {
    Path::new(&env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(1)
        .unwrap()
        .to_path_buf()
}

/// Check for compilation errors before building
fn check_compilation(sh: &Shell) -> Result<()> {
    println!("🔍 Checking for compilation errors...\n");

    let _dir = sh.push_dir(project_root());

    // Try to build the workspace
    let result = cmd!(sh, "cargo build --workspace --all-features")
        .ignore_status()
        .run();

    match result {
        Ok(_) => {
            println!("✅ No compilation errors found!");
            println!("Ready to build pycfdrs wheels.");
            Ok(())
        }
        Err(_) => {
            println!("❌ Compilation errors detected!");
            println!("\nCommon issues:");
            println!("  - Missing ToPrimitive trait bounds in cfd-1d/cfd-3d");
            println!("  - Type conversion issues in generic implementations");
            println!("\nFix these errors before building pycfdrs:");
            println!("  cargo build --workspace --all-features");
            bail!("Compilation failed")
        }
    }
}

/// Build pycfdrs wheel with maturin
fn build_wheel(sh: &Shell, release: bool) -> Result<()> {
    println!("🔨 Building pycfdrs wheel...\n");

    // First check compilation
    if let Err(_) = check_compilation(sh) {
        println!("\n⚠️  Fix compilation errors before building wheels.");
        println!("Run: cargo xtask check");
        return Ok(());
    }

    let pycfdrs_dir = project_root().join("crates").join("pycfdrs");
    let _dir = sh.push_dir(&pycfdrs_dir);

    // Check if maturin is installed
    if cmd!(sh, "maturin --version").ignore_status().run().is_err() {
        println!("📦 Installing maturin...");
        cmd!(sh, "pip install maturin").run()?;
    }

    // Build wheel
    if release {
        println!("Building release wheel...");
        cmd!(sh, "maturin build --release").run()?;
    } else {
        println!("Building debug wheel...");
        cmd!(sh, "maturin build").run()?;
    }

    // Find the wheel
    let wheel_dir = pycfdrs_dir.join("target").join("wheels");
    if wheel_dir.exists() {
        println!("\n✅ Wheel built successfully!");
        println!("Location: {}", wheel_dir.display());

        // List wheels
        if let Ok(entries) = std::fs::read_dir(&wheel_dir) {
            println!("\nAvailable wheels:");
            for entry in entries.flatten() {
                if entry.path().extension().and_then(|s| s.to_str()) == Some("whl") {
                    println!("  - {}", entry.file_name().to_string_lossy());
                }
            }
        }
    }

    Ok(())
}

/// Setup Python virtual environment
fn setup_venv(sh: &Shell, python: &str, venv_path: &Path) -> Result<()> {
    println!("🐍 Setting up Python virtual environment...\n");

    let _dir = sh.push_dir(project_root());

    // Check if venv already exists
    if venv_path.exists() {
        println!(
            "⚠️  Virtual environment already exists at: {}",
            venv_path.display()
        );
        println!("Remove it first to recreate: cargo xtask clean");
        return Ok(());
    }

    // Check Python version
    println!("Checking Python version...");
    let version_output = cmd!(sh, "{python} --version").read()?;
    println!("Using: {}", version_output);

    // Create virtual environment
    println!("\nCreating virtual environment...");
    cmd!(sh, "{python} -m venv {venv_path}").run()?;

    println!("✅ Virtual environment created at: {}", venv_path.display());
    println!("\nActivate with:");

    #[cfg(windows)]
    println!("  .venv\\Scripts\\activate");

    #[cfg(not(windows))]
    println!("  source .venv/bin/activate");

    Ok(())
}

/// Install Python dependencies
fn install_deps(sh: &Shell, venv_path: &Path, with_fenics: bool) -> Result<()> {
    println!("📦 Installing Python dependencies...\n");

    let _dir = sh.push_dir(project_root());

    // Get pip path in venv
    #[cfg(windows)]
    let pip = venv_path.join("Scripts").join("pip.exe");

    #[cfg(not(windows))]
    let pip = venv_path.join("bin").join("pip");

    if !pip.exists() {
        bail!("Virtual environment not found. Run: cargo xtask setup-venv");
    }

    let pip_str = pip.to_str().unwrap();

    // Upgrade pip
    println!("Upgrading pip...");
    cmd!(sh, "{pip_str} install --upgrade pip").run()?;

    // Install maturin
    println!("\nInstalling maturin...");
    cmd!(sh, "{pip_str} install maturin").run()?;

    // Install validation dependencies
    let requirements = project_root().join("validation").join("requirements.txt");
    if requirements.exists() {
        println!("\nInstalling validation dependencies...");
        cmd!(sh, "{pip_str} install -r {requirements}").run()?;
    }

    // Install pycfdrs wheel if available
    let wheel_dir = project_root()
        .join("crates")
        .join("pycfdrs")
        .join("target")
        .join("wheels");
    if wheel_dir.exists() {
        if let Ok(entries) = std::fs::read_dir(&wheel_dir) {
            for entry in entries.flatten() {
                if entry.path().extension().and_then(|s| s.to_str()) == Some("whl") {
                    println!("\nInstalling pycfdrs...");
                    let wheel = entry.path();
                    cmd!(sh, "{pip_str} install --force-reinstall {wheel}").run()?;
                    println!("✅ pycfdrs installed!");
                    break;
                }
            }
        }
    } else {
        println!("⚠️  No pycfdrs wheel found. Run: cargo xtask build-wheel");
    }

    if with_fenics {
        println!("\n⚠️  FEniCS installation requires conda.");
        println!("Run: cargo xtask install-fenics --use-conda");
    }

    println!("\n✅ Dependencies installed!");

    Ok(())
}

/// Install FEniCS for validation
fn install_fenics(sh: &Shell, use_conda: bool) -> Result<()> {
    println!("🧮 Installing FEniCS for validation...\n");

    if use_conda {
        // Check if conda is available
        if cmd!(sh, "conda --version").ignore_status().run().is_err() {
            println!("❌ Conda not found!");
            println!("\nFEniCS requires conda. Install from:");
            println!("  https://docs.conda.io/en/latest/miniconda.html");
            println!("\nOr use Docker:");
            println!("  docker run -ti quay.io/fenicsproject/stable:latest");
            return Ok(());
        }

        println!("Creating conda environment 'cfdrs-validation'...");
        println!("This may take several minutes...\n");

        // Create conda environment with FEniCS
        let result = cmd!(
            sh,
            "conda create -n cfdrs-validation -c conda-forge fenics python=3.11 -y"
        )
        .ignore_status()
        .run();

        match result {
            Ok(_) => {
                println!("\n✅ FEniCS installed in conda environment!");
                println!("\nActivate with:");
                println!("  conda activate cfdrs-validation");
                println!("\nThen install pycfdrs:");
                println!(
                    "  cargo xtask install-deps --venv $(conda info --base)/envs/cfdrs-validation"
                );
            }
            Err(_) => {
                println!("⚠️  Conda installation failed or environment already exists.");
                println!("To recreate:");
                println!("  conda env remove -n cfdrs-validation");
                println!("  cargo xtask install-fenics");
            }
        }
    } else {
        println!("⚠️  FEniCS pip installation is experimental and may not work.");
        println!("Conda installation is strongly recommended.");
        println!("\nRun with --use-conda flag:");
        println!("  cargo xtask install-fenics --use-conda");
    }

    Ok(())
}

/// Run validation suite
fn validate(sh: &Shell, venv_path: &Path, plot: bool, quick: bool, category: &str) -> Result<()> {
    println!("🧪 Running validation suite...\n");

    let _dir = sh.push_dir(project_root().join("validation"));

    // Get Python path in venv
    #[cfg(windows)]
    let python = venv_path.join("Scripts").join("python.exe");

    #[cfg(not(windows))]
    let python = venv_path.join("bin").join("python");

    if !python.exists() {
        bail!("Virtual environment not found. Run: cargo xtask setup-venv");
    }

    let python_str = python.to_str().unwrap();

    // Check if pycfdrs is installed
    let check_result = cmd!(sh, "{python_str} -c \"import pycfdrs\"")
        .ignore_status()
        .run();

    if check_result.is_err() {
        println!("❌ pycfdrs not installed in virtual environment!");
        println!("\nRun:");
        println!("  cargo xtask build-wheel");
        println!("  cargo xtask install-deps");
        return Ok(());
    }

    // Build command
    let mut args = vec!["run_all_validations.py"];

    if plot {
        args.push("--plot");
    }

    if quick {
        args.push("--quick");
    }

    args.push("--category");
    args.push(category);

    // Run validation
    println!("Running: python {}\n", args.join(" "));
    cmd!(sh, "{python_str} {args...}").run()?;

    println!("\n✅ Validation complete!");
    println!("Report: validation/reports/validation_summary.md");

    Ok(())
}

/// Complete workflow
fn run_all(sh: &Shell, with_fenics: bool, plot: bool) -> Result<()> {
    println!("🚀 Running complete validation workflow...\n");

    // 1. Check compilation
    println!("Step 1/6: Checking compilation...");
    if let Err(e) = check_compilation(sh) {
        println!("\n❌ Workflow stopped due to compilation errors.");
        println!("Fix errors and run: cargo xtask all");
        return Err(e);
    }
    println!();

    // 2. Build wheel
    println!("Step 2/6: Building pycfdrs wheel...");
    build_wheel(sh, true)?;
    println!();

    // 3. Setup venv
    println!("Step 3/6: Setting up virtual environment...");
    let venv_path = PathBuf::from(".venv");
    if !venv_path.exists() {
        setup_venv(sh, "python", &venv_path)?;
    } else {
        println!("Virtual environment already exists.");
    }
    println!();

    // 4. Install dependencies
    println!("Step 4/6: Installing dependencies...");
    install_deps(sh, &venv_path, with_fenics)?;
    println!();

    // 5. Install FEniCS (optional)
    if with_fenics {
        println!("Step 5/6: Installing FEniCS...");
        install_fenics(sh, true)?;
        println!();
    } else {
        println!("Step 5/6: Skipping FEniCS installation (use --with-fenics to enable)");
        println!();
    }

    // 6. Run validation
    println!("Step 6/6: Running validation suite...");
    validate(sh, &venv_path, plot, false, "analytical")?;

    println!("\n🎉 Complete workflow finished successfully!");
    println!("\nNext steps:");
    println!("  - Review report: validation/reports/validation_summary.md");
    if plot {
        println!("  - View plots: validation/reports/figures/");
    }
    println!("  - Run full validation: cargo xtask validate --plot");

    Ok(())
}

/// Clean build artifacts and virtual environment
fn clean(sh: &Shell) -> Result<()> {
    println!("🧹 Cleaning build artifacts...\n");

    let _dir = sh.push_dir(project_root());

    // Clean cargo
    println!("Cleaning Rust build artifacts...");
    cmd!(sh, "cargo clean").run()?;

    // Remove pycfdrs wheels
    let wheel_dir = project_root().join("crates").join("pycfdrs").join("target");
    if wheel_dir.exists() {
        println!("Removing pycfdrs wheels...");
        std::fs::remove_dir_all(&wheel_dir).context("Failed to remove wheel directory")?;
    }

    // Remove virtual environment
    let venv_path = project_root().join(".venv");
    if venv_path.exists() {
        println!("Removing virtual environment...");
        std::fs::remove_dir_all(&venv_path).context("Failed to remove virtual environment")?;
    }

    // Remove validation reports
    let reports_dir = project_root().join("validation").join("reports");
    if reports_dir.exists() {
        println!("Removing validation reports...");
        std::fs::remove_dir_all(&reports_dir).context("Failed to remove reports")?;
    }

    println!("\n✅ Clean complete!");

    Ok(())
}
