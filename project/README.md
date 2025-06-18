# Monte Carlo Electron Transport Simulation

This project simulates electron transport in gas mixtures using a Monte-Carlo approach to model collisions. It supports various gas species and computes transport and reaction properties based on cross-section data.

## Project Structure

- `src/` - Source code
  - `core/MonteCarlo.cpp, .hpp` - Main Monte Carlo class
  - `main/MC_singlerun.cpp` - main program to run a simulation
  - `utils/` - Utility classes for molar mass computation, cross sections, energy data, flux data, bulk data, and reaction rates
- `include/` - Header files
- `data` - Data files
  - `Xsec` - Cross-section data files for different species
  - `config` - File .json with simulation data
- `tests/` - Test codes

## Build

This project uses a Makefile for building. Available commands:

```bash
# Build the executable (default)
make

# Run the simulation using data in 'data/config/simulation.json'
make run

# Build with debug flags
make debug

# Remove build directory
make clean

# Show help with available targets
make help
```

The executable will be created in the `build/` directory as `mc_sim`.

## Usage

1. Build the executable with `make`.
2. Edit `data/config/simulation.json` to set simulation parameters and gas mixtures.
3. Run the simulation using `make run`.
4. Outputs will printed to the console roughly every 1e6 collisions.
5. Final results will be saved in a .txt file in folder `./results/`.
6. When done, clean the build directory with `make clean`.

**Note:** For multiple simulation runs, the executable does not require rebuilding. Simply repeat steps 2-5 for each simulation and run `make clean` only upon completion of all simulations.

## Extending

- Add new gas species by placing cross-section data files in `data/Xsec/`.
- Introduce a solver for the Poisson equation to allow the computation of a non-uniform, non-constant electric field.
- Improve the speed of randomic collisions simulations with a parallel implementation.

---