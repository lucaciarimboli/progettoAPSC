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

    // Carica configurazione
    ConfigParser config;
    
    std::string config_file = "data/config/simulation.json";
    if (argc > 1) {
        config_file = argv[1];  // use input file
    }
    
    if (!config.loadFromFile(config_file)) {
        std::cout << "Warning: using built-in default parameters" << std::endl;
    } else {
        std::cout << "Loaded simulation data from: " << config_file << std::endl;
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
    // flag to indicate the will of saving results in a file
    const bool save_in_file = config.save_results;

    std::cout << "\nStarting simulation...\n" << std::endl;

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
    if (duration >= 3600) {
        const int hours = duration / 3600;
        const int minutes = (duration % 3600) / 60;
        std::cout << "Simulation ended in " << hours << " hours " << minutes << " minutes\n";
    } else {
        const int minutes = duration / 60;
        std::cout << "Simulation ended in " << minutes << " minutes\n";
    }

     // Save results when simulation ends
    if(save_in_file){
        MC.saveResults(duration);
    }


    // FOR DEBUGGING PURPOSES:
    const std::vector<MeanData>& meandata = MC.get_mean_data();
    const std::vector<double>& time = MC.get_time_vector();

    std::ofstream file("tests/rz_data.csv");
    file << "time" << "r_z" << "var_z" << "\n";
    for(size_t i = 0; i < meandata.size(); i++) {
        file << time[i] << "," << meandata[i].get_position()[2] << "," << meandata[i].get_variance()[2] << "\n";
    }
    file.close();
    
    return 0;
}