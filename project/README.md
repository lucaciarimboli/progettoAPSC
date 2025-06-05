# Monte Carlo Electron Transport Simulation

This project simulates electron transport in gas mixtures using a Monte-Carlo approach to model collisions. It supports various gas species and computes transport and reaction properties based on cross-section data.

## Project Structure

- `src/` - Source code
  - `core/MonteCarlo.cpp, .hpp` - Main Monte Carlo class
  - `main/MC_singlerun.cpp` - main program to run a simulation
  - `utils/` - Utility classes for cross-sections, molar mass, energy, flux, bulk data, and reaction rates
- `include/` - Common header
- `data/Xsec/` - Cross-section data files for different species
- `tests/` - Test codes

## Build

To build the project, use your preferred C++ build system (e.g., CMake or Make). Make sure to include all source files in `src/` and set the include path to `include/`.

## Usage

1. Edit `src/main/MC_singlerun.cpp` to set simulation parameters and gas mixtures.
2. Build and run the executable.
3. Output and results will be printed to the console.

## Extending

- Add new gas species by placing cross-section data files in `data/Xsec/`.
- Introduce a solver for the Poisson equation to allow the computation of a non-uniform, non-constant electric field.
- Improve the speed of randomic collisions simulations with a parallel implementation.

---