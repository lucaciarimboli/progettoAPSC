# Monte Carlo Electron Transport Simulation

This project simulates electron transport in a low temperature plasma using a Monte-Carlo, null-collision approach to model
electron-neutral collisions with the backgorund gas particles.
The code supports various gas species and computes transport parameters, reaction rates and energy distribution of the electrons
based on cross-section data (in the LXCat format).

## Project Structure

- `src/` - Source code
  - `core/MonteCarlo.cpp, .hpp` - Main Monte Carlo class
  - `main/MC_singlerun.cpp` - main program to run a simulation
  - `utils/` - Utility classes for molar mass computation, cross sections, energy data, flux data, bulk data, and reaction rates
- `include/` - Header files
- `data` - Data files
  - `Xsec` - Cross-section data files for different species
  - `config` - File .json with simulation data
- `tests/` - Some test codes for implementation

## Build

This project uses a Makefile for building. Available commands are:

```bash
# Build the executable (default)
make

# Run the simulation using data in 'data/config/simulation.json'
make run

# Generate Doxygen documentation
make docs

# Remove build and docs directories
make clean

# Build and run with profiler enabled
make profile

# Show help with available targets
make help
```

The executable will be created in the `build/` directory as `mc_sim`.
The Doxygen documentation files will be created in the `docs/` directory.
In case of gprof profiler enabled, the output file is saved in the `profiler/` directory as `report.txt`

## Usage

1. Build the executable with `make`.
2. Edit `data/config/simulation.json` to set gas mixture and simulation parameters.
3. Run the simulation using `make run`.
4. Partial data will be printed to the console roughly every 1e6 collisions.
5. Final results will be saved in a .txt file in `results/`.
6. Repeat steps 2-5 for running multiple simulations

## Extending

- Add new gas species by including cross-section data files in `data/Xsec/`.
- Implement a Poisson equation solver to allow simulations with a non-uniform, non-constant electric field.

---