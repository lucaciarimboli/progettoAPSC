#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <chrono>
#include "core/MonteCarlo.hpp"
#include "ConfigParser.hpp"

int main(int argc, char* argv[]) {

    // Start the clock
    auto start = std::chrono::high_resolution_clock::now();
    
    //-----------------------------------------//
    //       IMPORT SIMULATION PARAMETERS:     //
    //-----------------------------------------//

    ConfigParser config;
    
    std::string config_file = "data/config/simulation.json";
    if (argc > 1) {
        config_file = argv[1];
    }
    
    if (config.load_from_file(config_file)) {
        std::cout << "\n Loaded simulation data from: " << config_file << std::endl;
    } else {
        std::cout << "\n Failed to load simulation data from: " << config_file << std::endl;
    }

    //------------------------------------------//
    //       INITIALIZE MonteCarlo OBJECT       //
    //------------------------------------------//

    MonteCarlo MC(
        config.gas_species, 
        config.gas_mixture, 
        config.EN_field, 
        config.pressure, 
        config.temperature,
        config.initial_electrons, 
        config.max_electrons, 
        config.energy_sharing,
        config.max_energy, 
        config.energy_step,
        config.drift_velocity_error, 
        config.diffusion_error,
        config.min_collisions, 
        config.max_collisions,
        config.initial_position, 
        config.initial_broadening,
        config.conserve_electrons, 
        config.isotropic_scattering
    );

    //------------------------------------------//
    //         RUN MONTE CARLO SIMULATION       //
    //------------------------------------------//
    
    // checks end of simulation: End =1 stops the simulation
    bool End = false;

    std::cout << " Starting simulation..." << std::endl;

    while( !End ) {
        // Perform a flight for all electrons without a collision:
        MC.freeFlight();
        // Calculate the collective data of the electron swarm:
        MC.collectMeanData();

        // If there is enough data at steady state, update accordingly the variables of interest:
        if(MC.get_count_sst() > 10) {
            MC.updateEnergyData();
            MC.updateFluxData();
            MC.updateBulkData();
            MC.updateReactionRates();
        }

        // Perform collisions:
        MC.updateCollisionMatrix();
        MC.performCollisions();

        // Check if the simulation has reached steady state:
        MC.checkSteadyState();

        // Print the meaningful data at the current time step:
        MC.printOnScreen();
        // Check if the simulation has converged:
        End = MC.endSimulation();
    }

    // Print the simulation time:
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

    const int minutes = duration / 60;
    const int seconds = duration % 60;
    std::cout << " Simulation ended in " << minutes << " minutes, " << seconds << " seconds\n";    

    // Save results when simulation ends
    MC.saveResults(duration);

    return 0;
}