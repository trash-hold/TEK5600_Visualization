# TEK5600 Visualization Exercises

This project automatically builds executables for each subdirectory in `exercises/` that contains a `.cpp` file.

## Current Exercises

### Field Lines (`field_lines/`)

- **`field_lines_benchmark_generators`**: Benchmarks different seed generation methods for field lines
- **`field_lines_display_line`**: Displays individual field lines
- **`field_lines_even_fieldlines_heatmap`**: Evenly spaced field lines with velocity magnitude heatmap (inferno colormap)
- **`field_lines_method_comp`**: Comparison of field line generation methods

### Integrator Comparison (`integrator_comp/`)

- **`integrator_comp_benchmark_integrators`**: Benchmarks Euler vs RK4 integrators
- **`integrator_comp_comp_isabel`**: Integrator comparison on Isabel dataset
- **`integrator_comp_comp_metsim`**: Integrator comparison on Metsim dataset

### LIC (Line Integral Convolution) (`lic/`)

- **`lic_lic_blending`**: LIC with velocity magnitude heatmap and multiple blending modes (Multiply, Screen, Soft Light)
- **`lic_lic_integrator_comparison`**: LIC performance comparison between Euler and RK4 integrators
- **`lic_simple_lic`**: Basic LIC visualization

## Directory Structure

```
exercises/
├── field_lines/
│   ├── benchmark_generators.cpp
│   ├── display_line.cpp
│   ├── even_fieldlines_heatmap.cpp
│   ├── method_comp.cpp
│   └── img/
├── integrator_comp/
│   ├── benchmark_integrators.cpp
│   ├── comp_isabel.cpp
│   ├── comp_metsim.cpp
│   └── img/
├── lic/
│   ├── lic_blending.cpp
│   ├── lic_integrator_comparison.cpp
│   ├── simple_lic.cpp
│   └── img/
└── README.md
```

## Building

From the build directory:
```bash
cmake ..
make
```

CMake automatically detects all `.cpp` files in `exercises/` subdirectories and builds executables named `{dirname}_{filename}`.

## Running

### Field Lines
```bash
./field_lines_benchmark_generators    # Benchmark seed generators
./field_lines_display_line            # Display field lines
./field_lines_even_fieldlines_heatmap # Field lines with heatmap
./field_lines_method_comp             # Method comparison
```

### Integrator Comparison
```bash
./integrator_comp_benchmark_integrators # Benchmark integrators
./integrator_comp_comp_isabel          # Compare on Isabel data
./integrator_comp_comp_metsim          # Compare on Metsim data
```

### LIC
```bash
./lic_simple_lic                 # Basic LIC
./lic_lic_integrator_comparison  # LIC integrator comparison
./lic_lic_blending               # LIC with blending modes (press 1/2/3 to switch)
```

## Adding New Exercises

1. Create a new subdirectory in `exercises/` (e.g., `exercises/new_exercise/`)
2. Add a `.cpp` file in that directory
3. Include necessary headers from `inc/`
4. Rebuild with `make` - the executable will be automatically created as `{dirname}_{filename}`

No manual CMakeLists.txt editing required! The builds are kept in the main `build/` folder.