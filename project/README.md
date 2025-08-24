# Monte Carlo Electron Transport Simulation

This project simulates electron transport in a low temperature plasma using a Monte-Carlo, null-collision approach to model
electron-neutral collisions with the backgorund gas particles.
The code supports various gas species and computes transport parameters, reaction rates and energy distribution of the electrons
based on cross-section data (in the LXCat format).

## Project Structure
```
.
|-- Doxyfile
|-- Makefile
|-- data
|   |-- Xsec
|   `-- config
|       `-- simulation.json
|-- include
|   |-- Common.hpp
|   |-- ConfigParser.hpp
|   |-- core
|   |   `-- MonteCarlo.hpp
|   |-- nlohmann
|   `-- utils
|       |-- BulkData.hpp
|       |-- CollisionData.hpp
|       |-- CrossSectionsData.hpp
|       |-- EnergyData.hpp
|       |-- FluxData.hpp
|       |-- MeanData.hpp
|       |-- MolMass.hpp
|       `-- ReactionRates.hpp
|-- results
`-- src
    |-- core
    |   `-- MonteCarlo.cpp
    |-- main
    |   `-- MC_singlerun.cpp
    `-- utils
        |-- BulkData.cpp
        |-- CollisionData.cpp
        |-- CrossSectionsData.cpp
        |-- EnergyData.cpp
        |-- FluxData.cpp
        |-- MeanData.cpp
        |-- MolMass.cpp
        `-- ReactionRates.cpp
```

## Build

This project uses a Makefile for building. Available targets are:

```bash
# Build the executable
make

# Run the simulation using data in 'data/config/simulation.json'
make run

# Run multiple simulations using data in 'data/config/simulation.json'
make multi_run

# Generate Doxygen documentation
make docs

# Remove build and docs directories
make clean

# Build and run with profiler enabled
make profile

# Show help with available targets
make help
```

The executable will be created in the `build/` directory with the name `mc_sim`.
The Doxygen documentation files will be created in the `docs/` directory.
In case of run with profiler enabled, the gprof report is saved in the `profiler/` directory as `report.txt`.

## How to perform a simulation

1. Build the executable with `make`.
2. Edit `data/config/simulation.json` to set gas mixture and simulation parameters.
3. Run the simulation using `make run`.
4. Partial data will be printed to the console roughly every 1e6 collisions.
5. Final results of the simulation will be saved in a .txt file in `results/`.

To run multiple simulations with the same data, after steps *1.* and *2.* call the target `make multi_run`.
After calling the target, the desired number of simulations should be typed in the shell as required.
The results of subsequent simulations will be saved in `results/` with progressively increasing naming, e.g. `results/result1.txt`, `results/result2.txt`, and so forth...

## How to set up data

The simulation data should be defined in the file `data/config/simulation.json`.
Before running a simulation, this file should be opened and the input data for the simulation should
be customized at will. Default typical data to perform a simulation are already inserted.
The required data in `simulation.json` are divided in:

- *Simulation data*:
  - `gas_species`: Gas species involved in the mixture. Every gas species is identified by its chemical formula, e.g. `["N2", "O2"]`.
  Currently supported species are: H2, H2O, N, N2, O, O2. To use other species, please download the cross sections data files from the [LXCat database](https://us.lxcat.net/data/set_type.php).
  - `gas_mixture`: fractions mixture components specified in `gas_species`. The number of fractions should be the same as the number of species.
  The fractions should sum to 1, if they do no, the code automatically modifies the last entry to make the total fraction equal 1.
  - `EN_field`: electric field intensity per electron, in Td.
  - `pressure`: gas pressure, in Pa.
  - `temperature`: gas temperature, in K.
  - `initial_electrons`: inital number of the electrons in the simulation.
  - `energy_sharing`: energy sharing factor for ionization events. It should be a value in between 0 and 1. 
  - `initial_position`: vector of initial center of mass of initial gaussian distributed electrons in x,y and z direction, e.g. `[0.0, 0.0, 0.0]`.
  - `initial_broadening`: vector of initial broadening of initial gaussian distributed electrons in x,y and z direction, e.g. `[0.0, 0.0, 0.0]`.

- *Energy grid data*:
  - `max_energy`: maximum electron energy, in eV.
  - `energy_step`: step for the uniform energy grid, in eV.

- *Tolerances*:
  - `drift_velocity_error`: tolerance threshold for the relative error in bulk velocity.
  - `diffusion_error`: tolerance threshold for the relative error in bulk diffusion coefficient.
  - `min_collisions`: minimum number of collisions to check for steady state condition.
  - `max_collisions`: maximum number of collisions allowed.
  - `max_electrons`: maximum electon population size allowed.

- *Flags*:
  - `conserve_electrons`: if true, enforces constant size of the electron population.
  - `isotropic_scattering`: if true, enforces isotropic scattering. 


## To do

- Add new gas species by including cross-section data files in `data/Xsec/`.
- Implement a Poisson equation solver to allow simulations with a non-uniform, non-constant electric field.

---
