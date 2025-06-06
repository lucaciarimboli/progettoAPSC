#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <chrono>
#include "core/MonteCarlo.hpp"

int main() {

    // Start the clock
    auto start = std::chrono::high_resolution_clock::now();
    
    //-----------------------------------------//
    //          SIMULATION PARAMETERS:         //
    //-----------------------------------------//

    // Gas species and their fractions: ( allowed species: "H2", "H2O", "N", "N2", "O", "O2" )
    //std::vector<std::string> gas = {"N2", "N", "O2", "O", "H2O", "H2"};
    //std::vector<double> mix = {0.8, 0.0, 0.2, 0.0, 0.0, 0.0};

    const std::vector<std::string> gas = {"N2", "O2"};
    const std::vector<double> mix = {0.8, 0.2};

    // E/N (in Td):
    const double EN = 5e3;
    // Pressure (in Pa):
    const double p = 101325;
    // Temperature (in K):
    const double T = 300;
    // Initial number of electrons:
    const unsigned long N0 = 1e4;
    // Maximum allowed number of electrons:
    const unsigned long Ne_max = 1e6;
    // Energy sharing factor for ionization (in interval [0,1]):
    const double W = 0.5;

    // Define the maximum energy level and the step for electrons energy grid (in eV):
    const double E_max = 1e5;
    const double dE = 0.1;   // then energy grid will be [0,E_max] with step "E_step"

    // Error tolerances:
    const double w_err = 0.01;    // for drift velocity  
    const double DN_err = 0.01;   // for diffusion constant

    // Minimum number of collisions for steady state:
    const unsigned long col_equ = 1e7;
    // Maximum number of collisions:
    const unsigned long col_max = 1e8;

    // Initial mean position of electrons (in m):
    const std::array<double, 3> pos_xyz = {0.0, 0.0, 0.0};
    // Initial broadening of electrons (in m):
    const std::array<double, 3> sigma_xyz = {0.0, 0.0, 0.0};

    // Flag for conserving electron number after ionization/attachment:
    const bool conserve = true;
    // Flag for isotropic scattering:
    const bool isotropic = true;

    //------------------------------------------//
    //       INITIALIZE MonteCarlo OBJECT       //
    //------------------------------------------//

    MonteCarlo MC
    (
        gas, mix, EN, p, T, N0, Ne_max, W, E_max, dE, w_err, DN_err,
        col_equ, col_max, pos_xyz, sigma_xyz, conserve, isotropic
    );
    

    //------------------------------------------//
    //         RUN MONTE CARLO SIMULATION       //
    //------------------------------------------//
    
    // checks end of simulation: End =1 stops the simulation
    bool End = false;


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
        MC.performCollision("ELASTIC");
        MC.performCollision("EXCITATION");
        MC.performCollision("IONIZATION");
        MC.performCollision("ATTACHMENT");

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
        std::cout << "Reaction rates updated in " << hours << " hours " << minutes << " minutes\n";
    } else {
        const int minutes = duration / 60;
        std::cout << "Reaction rates updated in " << minutes << " minutes\n";
    }
    
    return 0;
}