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

    // Constructor:
    MonteCarlo( const std::vector<std::string> & gas, const std::vector<double> & mix,
                const double EN, const double p, const double T, const unsigned N0, const unsigned Ne_max,
                const double W, const double E_max, const double dE,
                const double w_err, const double DN_err, const unsigned col_equ, const unsigned col_max,
                const std::array<double,3> & pos_xyz, const std::array<double,3> & sigma_xyz,
                const bool conserve, const bool isotropic);
    
    // Public Methods:
    void freeFlight();
    void collectMeanData();
    void updateEnergyData();
    void updateFluxData();
    void updateBulkData();
    void updateReactionRates();
    void updateCollisionMatrix();
    void performCollisions();
    void checkSteadyState();
    bool endSimulation();
    void printOnScreen();
    void saveResults(const int64_t duration) const;

    // Getters:
    const unsigned int get_count_sst() const { return count_sst; }

    //---------------------------------------------------------------------//
    //------------------- FOR DEBUGGING PURPOSES --------------------------//
    const std::vector<MeanData>& get_mean_data() const { return mean; }
    const std::vector<double>& get_time_vector() const { return t; }
    //---------------------------------------------------------------------//
    //---------------------------------------------------------------------//

private:

    // Class Members:
    const unsigned  N0;     // initial electron population   
    const double N;         // gas density in m^-3
    
    const std::vector<std::string> gas; // subformulas of gas species
    std::vector<double> mgas;           // mass of gas species (in kg)
    std::vector<double> mix;            // fractions of individual species
            
    CrossSectionsData Xsec;         // cross-section data for gas species
                    
    const double w_err;             // error tolerance of mean energy (default: 1%)
    const double DN_err;            // error tolerance of diffusion constant (default: 1%)
    const unsigned Ne_max;          // maximum allowed number of electrons
    const unsigned long col_equ;    // minimum # of collisions for equilibrium (default: 20e6)
    const unsigned long col_max;    // max # of collisions (default: 20e6) 
    const bool conserve;            // force electron conservation (true), allow gain/loss (false)
    const bool isotropic;           // isotropic scattering (true) or not (false)

    const double W;                 // energy sharing in ionizing collision
    const double sqrtW;             // std::sqrt(W) -> pre computation
    const double sqrt1_W;           // std::sqrt(1-W) -> pre computation
    double E_max;                   // maximum electron energy
    const double EN;                // Electric field per-electron [Td]

    double T_sst = 0.0;             // equilibrium time
    unsigned int count_sst = 0;     // counts the number of steady state time steps
    int line = 1;   // line number in output file:

    std::vector<double> t;                          // vector of times [s]
    double dt;                                      // current time step [s]
    double t_total = 0.0;                           // sum of all times for all electrons [m]
    mc::MATRIX r;                                   // 3D positions of electrons [m]
    mc::MATRIX v;                                   // current velocity of electrons [m/s]
    std::array<double,3> a;                         // acceleration of electrons [m/s^2] (const & uniform)
    mc::MATRIX v_int;                               // time-integrated velocity
    mc::MATRIX v2_int;                              // time-integrated squared velocity
    mc::MATRIX r_cations;                           // 3D positions of cations; [m]
    mc::MATRIX r_anions;                            // 3D positions of anions;  [m]

    CollisionData C;                                // Collision indeces, Mass and Loss vectors
    unsigned long collisions = 0;                   // total number collisions occourred

    std::vector<MeanData> mean;     // temporal mean data of electron swarm
    EnergyData E;                   // energy data
    BulkData bulk;                  // bulk transport data
    FluxData flux;                  // flux transport data
    RateDataConv rates_conv;        // reaction rates computed by convolution
    RateDataCount rates_count;      // reaction rates computed by counting

    unsigned int converge = 0;      // Status
            
    double EnergyLossElastic    = 0;    // Heat energy loss in elastic collisions
    double EnergyLossInelastic  = 0;    // Heat energy loss in inelastic collisions
    double EnergyLossIonization = 0;    // heat energy loss in ionization collisions

    std::default_random_engine gen;             // Random number generator
    std::uniform_real_distribution<> randu;     // Uniform distribution for random numbers
    std::normal_distribution<double> randn;     // Normal distribution for random numbers


    // Private Methods:
    void checkFractionSum();
    void mass_in_kg();
    void initialParticles(const std::array<double,3> & pos_xyz, const std::array<double,3> & sigma_xyz);
    std::pair<double,double> velocity2energy(const std::array<double,3> & v) const;
    std::array<double, 3> cross_product(const std::array<double, 3>& a, const std::array<double, 3>& b) const;
    void elasticCollision(const std::vector<size_t> & ind, const std::vector<double> & Mass);
    void inelasticCollision(const std::vector<size_t> & ind, const std::vector<double> & Loss);
    void ionizationCollision(const std::vector<size_t> & ind, const std::vector<double> & Loss);
    void attachmentCollision(const std::vector<size_t> & ind);
};

#endif  // MONTECARLO_HPP