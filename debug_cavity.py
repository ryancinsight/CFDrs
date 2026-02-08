
import pycfdrs
import numpy as np

def debug_cavity():
    print("Initializing solver (33x33)...")
    # Low resolution to avoid memory/stability issues during debug
    solver = pycfdrs.PyCavitySolver2D(nx=33, ny=33, cavity_size=1.0, lid_velocity=1.0, reynolds=100.0)
    
    print("Running solver...")
    try:
        result = solver.solve()
        print(f"Solver finished. L2 Error: {result.l2_error:.4f}")
    except Exception as e:
        print(f"Solver failed: {e}")

if __name__ == "__main__":
    debug_cavity()
