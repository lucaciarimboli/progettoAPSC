#include <iostream>
#include <vector>
#include <array>
#include <string>
#include "MonteCarlo.hpp"

int main() {
    
    //-----------------------------------------//
    //          SIMULATION PARAMETERS:         //
    //-----------------------------------------//

    // Gas species and their fractions: ( allowed species: "H2", "H2O", "N", "N2", "O", "O2" )
    std::vector<std::string> gas = {"N2", "N", "O2", "O", "H2O", "H2"};
    std::vector<double> mix = {0.8, 0.0, 0.2, 0.0, 0.0, 0.0};

    // E/N (in Td):
    double EN = 5e3;
    // Pressure (in Pa):
    double p = 101325;
    // Temperature (in K):
    double T = 300;
    // Initial number of electrons:
    unsigned N0 = 10000;
    // Maximum allowed number of electrons:
    unsigned Ne_max = 1e6;
    // Energy sharing factor for ionization (in interval [0,1]):
    double W = 0.5;

    // Define the maximum energy level and the stp for electrons energy grid (in eV):
    double E_max = 100;
    double dE = 1;   // then energy grid will be [0,E_max] with step "E_step"

    // Error tolerances:
    double w_err = 0.01;    // for drift velocity  
    double DN_err = 0.01;   // for diffusion constant

    // Minimum number of collisions for steady state:
    unsigned col_equ = 1e7;
    // Maximum number of collisions:
    unsigned col_max = 1e8;

    // Initial mean position of electrons (in m):
    std::array<double, 3> pos_xyz = {0.0, 0.0, 0.0};
    // Initial broadening of electrons (in m):
    std::array<double, 3> sigma_xyz = {0.0, 0.0, 0.0};

    // Flag for conserving electron number after ionization/attachment:
    bool conserve = true;
    // Flag for isotropic scattering:
    bool isotropic = true;

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
        if(MC.countSteadyState()){
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
    
    return 0;
}