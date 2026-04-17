# TEK5600 Visualization

A C++ project for scientific data visualization techniques, focusing on vector field visualization methods including field lines and Line Integral Convolution (LIC).

## Features

- **Field Lines**: Implementation of Euler and RK4 integration methods with various seeding strategies
- **LIC (Line Integral Convolution)**: Texture-based flow visualization with multiple blending modes
- **Interactive Visualizations**: Real-time parameter adjustment and comparison tools
- **Scientific Data Support**: HDF5 file format for vector field datasets

## Dependencies

- C++17 compiler
- CMake 3.15+
- SFML 3.0 (graphics library)
- HDF5 (data format library)

## Building

```bash
mkdir build
cd build
cmake ..
make
```

## Project Structure

```
├── CMakeLists.txt          # Build configuration
├── inc/                    # Header files (primitives, algorithms)
├── src/                    # Implementation files
├── exercises/              # Executable implementations
│   ├── field_lines/        # Field line visualizations
│   ├── integrator_comp/    # Integration method comparisons
│   └── lic/                # LIC implementations
├── data/                   # Sample datasets (HDF5 format)
└── build/                  # Build output (generated)
```

## Running

After building, executables are created in the `build/` directory. See `exercises/README.md` for detailed descriptions of each visualization tool.

Example:
```bash
cd build
./lic_lic_blending  # Interactive LIC with blending modes
```

## Datasets

The project includes two sample vector field datasets:
- `isabel_2d.h5`: Atmospheric flow data
- `metsim1_2d.h5`: Meteorological simulation data