# Exercise Executables

This project automatically builds executables for each subdirectory in `exercises/` that contains a `main.cpp` file.

## Current Exercises

- **`lic_exercise`**: Runs Line Integral Convolution (LIC) visualization on the metsim data file.
- **`field_lines_exercise`**: Generates and displays field lines using even seeding on the isabel data file.
- **`ex1_test`**: Original combined executable (LIC visualization).

## Directory Structure

```
exercises/
├── lic/
│   └── main.cpp          # LIC visualization
├── field_lines/
│   └── main.cpp          # Field lines visualization
└── README.md
```

## Building

From the build directory:
```bash
cmake ..
make
```

CMake automatically detects all subdirectories in `exercises/` and builds executables named `{dirname}_exercise`.

## Running

```bash
./lic_exercise          # LIC visualization
./field_lines_exercise  # Field lines visualization
./ex1_test             # Original LIC
```

## Adding New Exercises

1. Create a new subdirectory in `exercises/` (e.g., `exercises/new_exercise/`)
2. Add a `main.cpp` file in that directory
3. Include necessary headers from `inc/`
4. Rebuild with `make` - the executable will be automatically created as `new_exercise_exercise`

No manual CMakeLists.txt editing required! The builds are kept in the main `build/` folder.