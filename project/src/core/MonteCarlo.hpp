#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <numeric>
#include <cmath>
#include <utility>
#include <random>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "utils/CrossSectionsData.hpp"
#include "utils/MolMass.hpp"
#include "utils/MeanData.hpp"
#include "utils/EnergyData.hpp"
#include "utils/FluxData.hpp"
#include "utils/BulkData.hpp"
#include "utils/ReactionRates.hpp"
#include "utils/CollisionData.hpp"

#include "Common.hpp"


class MonteCarlo
{

public:
    // Constructor
    MonteCarlo( const std::vector<std::string> & gas, const std::vector<double> & mix,
                const double EN, const double p, const double T, const unsigned N0, const unsigned Ne_max,
                const double W, const double E_max, const double dE,
                const double w_err, const double DN_err, const unsigned col_equ, const unsigned col_max,
                const std::array<double,3> & pos_xyz, const std::array<double,3> & sigma_xyz,
                const bool conserve, const bool isotropic);
    
    // Performs non-collissional flight for electrons in electric field
    void freeFlight();
    // Gets electron collective data: mean position, mean
    // broadening in x,y and-direction, mean kinetic energy, electron number
    // and total electron current
    void collectMeanData();
    // Calculates mean energy and EEDF data after steady state was reached
    void updateEnergyData();
    // Calculates flux data after steady state was reached
    void updateFluxData();
    // Calculates bulk data after steady state was reached
    void updateBulkData();
    // Calculates reaction rates after steady state was reached
    void updateReactionRates();
    // Decides which collision will happen for each electron
    void updateCollisionMatrix();
    // Perform collisions:
    void performCollision(const std::string & type);
    // Checks if the simulation has reached steady state.
    // If so, it updates the time step and the number of collisions
    void checkSteadyState();
    // Stops the simulation
    bool endSimulation();
    // Prints the results of the simulation every iteration
    void printOnScreen();
    // Prints the final results of the simulation to a file
    void saveResults(const int64_t duration) const;

    // Getters:
    const unsigned int get_count_sst() const { return count_sst; }


    //------------------- FOR DEBUGGING PURPOSES --------------------------//
    const std::vector<MeanData>& get_mean_data() const { return mean; }
    const std::vector<double>& get_time_vector() const { return t; }
    //---------------------------------------------------------------------//

private:

    // number of initial electrons used in MC calculation
    const unsigned  N0;        
    // gas density in m^-3
    const double N;
    
    // cell array of subformula of gas species
    const std::vector<std::string> gas;
    // cell array of mass of gas species (in kg)
    std::vector<double> mgas;
    // fractions of individual species in the gas as a vector
    std::vector<double> mix;
            
    // cross sections data
    CrossSectionsData Xsec;
                    
    // tolerance in error of drift velcocity (default: 1%)
    const double w_err;
    // tolerance in error of diffusion constant (default: 1%)
    const double DN_err;
    // maximum allowed number of electorns
    const unsigned Ne_max;
    // number of collisions until equilibrium (default: 20e6)
    const unsigned long col_equ;
    // number of collisions at which simulation ends (default: 20e6) 
    const unsigned long col_max;
    // conserve (1) electron number after ionizatzion/attachment or not (0)
    const bool conserve;    
    // (1) isotropic, (0) non-isotropic scattering according to Vahedi et al.
    const bool isotropic;
    // energy sharing in ionizing collision
    const double W;
    // maximum electron energy
    double E_max;
    // collision counter
    // unsigned int counter = 0;
    // equilibrium time
    double T_sst = 0.0;
    // counts the number of steady state time steps
    unsigned int count_sst = 0;

    // line number in output file:
    int line = 1;
    // computation time:
    // elapsedTime;

    // current time [s]
    std::vector<double> t;
    // current time step dt
    double dt;
    // sum of all times for all electrons:
    double t_total = 0.0;
    // current position of electrons, cations and anions (order important)
    // r has 3 vectors: one for each type of particles.
    // e.g. vector for e contains the position (array<double,3>) of every electron;
    // x-coordinate of i-esim electron can be accessed through: r[ELECTRONS][i][0].
    std::array<mc::MATRIX,mc::PARTICLES_TYPES> r;
    // current velocity of electrons
    mc::MATRIX v;
    // acceleration of electrons (constant & uniform as the E field)
    std::array<double,3> a;
    // current time-integrated velocity
    mc::MATRIX v_int;
    // current time-integrated velocity-squared
    mc::MATRIX v2_int;
    // Collision matrix and collision indeces:
    CollisionData C;
    // total number of all real collisions that happend
    unsigned long collisions = 0;

    // temporal mean data of electron swarm
    std::vector<MeanData> mean; 
    // energy data
    EnergyData E;
    // bulk transport data
    BulkData bulk;
    // flux transport data
    FluxData flux;
    // reaction rates
    RateDataConv rates_conv;
    RateDataCount rates_count;
            
    // electric field function in x-direction
    double E_x;
    // electric field function in y-direction
    double E_y;
    // electric field function in z-direction
    double E_z;
            
    //Heat energy partition
    double EnergyLossElastic    = 0;
    double EnergyLossInelastic  = 0;
    double EnergyLossIonization = 0;
            
    // Status
    unsigned int converge = 0;

    // Random number generator:
    std::default_random_engine gen;
    // Uniform distribution for random numbers:
    std::uniform_real_distribution<> randu;


    // PRIVATE METHODS:

    // Checks if sum of gas fractions is equal to 1
    // if not the case: last entry of mix will be corrected
    void checkFractionSum();
    // Computes masses [kg] of the gas species
    void mass_in_kg();
    // Sets initial position and velocity of electrons
    void initialParticles(const std::array<double,3> & pos_xyz, const std::array<double,3> & sigma_xyz);
    // Sets E parallel to z-direction --> replaces "solvePoissin_3D()"
    void set_E(const double EN);
    // Calculates absolute value of velocity and energy in eV
    std::pair<double,double> velocity2energy(const std::array<double,3> & v) const;
    // Computes cross product of two 3D vectors (needed in methods that perform collisions):
    std::array<double, 3> cross_product(const std::array<double, 3>& a, const std::array<double, 3>& b) const;
    // Performs elastic collision (isotropic or non-isotropic)
    void elasticCollision(const std::vector<size_t> & ind, const std::vector<double> & Mass);
    // Performs inelastic collision (isotropic or non-isotropic)
    void inelasticCollision(const std::vector<size_t> & ind, const std::vector<double> & Loss);
    // Performs ionization collision for electrons
    void ionizationCollision(const std::vector<size_t> & ind, const std::vector<double> & Loss);
    // Performs attachment collision for electrons
    void attachmentCollision(const std::vector<size_t> & ind);
};

#endif  // MONTECARLO_HPP