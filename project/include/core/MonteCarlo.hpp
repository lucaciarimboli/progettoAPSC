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


/**
 * @class MonteCarlo
 * @brief Simulates electron transport in LTPs using the null-collision Monte Carlo method.
 *
 * This class models the dynamics of electrons in a low-temperature, non-Maxwellian plasma subjected to an external electric field.
 * It uses a Monte-Carlo approach with the null-collision technique to simulate interactions 
 * between electrons and background gas molecules.
 * Under steady-state conditions, it allows the computation and update of transport coefficients, 
 * reaction rates and energy distributions.
 * 
 * The MonteCarlo class is designed to be self-contained: it autonomously handles
 * the update of its members by calling their methods or extracting values from them.
 * As a result, the main file only needs to define a MonteCarlo object, 
 * without directly interacting with internal members data structures.
 * 
 */
class MonteCarlo
{

public:

    /**
     * @brief Constructor.
     * @param gas List of gas species.
     * @param mix Corresponding fractions of each gas species.
     * @param EN Reduced electric field [Td].
     * @param p Pressure [Pa].
     * @param T Temperature [K].
     * @param N0 Initial number of electrons.
     * @param Ne_max Maximum allowed number of electrons.
     * @param W Energy sharing factor for ionizing collisions.
     * @param E_max Maximum electron energy [eV].
     * @param dE Uniform step for energy grid [eV].
     * @param w_err Tolerance on bulk velocity.
     * @param DN_err Tolerance on diffusion coefficient.
     * @param col_equ Minimum number of collisions for equilibrium.
     * @param col_max Maximum number of total collisions allowed.
     * @param pos_xyz Mean initial position of electrons [m].
     * @param sigma_xyz Standard deviation of initial electron position [m].
     * @param conserve Flag for enforcing conservation of electrons.
     * @param isotropic Flag for isotropic scattering.
     */
    MonteCarlo( const std::vector<std::string> & gas, const std::vector<double> & mix,
                const double EN, const double p, const double T, const unsigned N0, const unsigned Ne_max,
                const double W, const double E_max, const double dE,
                const double w_err, const double DN_err, const unsigned col_equ, const unsigned col_max,
                const std::array<double,3> & pos_xyz, const std::array<double,3> & sigma_xyz,
                const bool conserve, const bool isotropic);
    
    // Public Methods:

    /**
     * @brief Simulates free flight of electrons between collisions.
     */
    void freeFlight();

    /**
     * @brief Collects center of mass, mean velocity and mean kinetic energy of the swarm.
     */    
    void collectMeanData();

    /**
     * @brief Updates energy distribution data at steady-state.
     */
    void updateEnergyData();

    /**
     * @brief Updates flux transport parameters at steady state.
     */
    void updateFluxData();

    /**
     * @brief Updates bulk transport parameters at steady state.
     */
    void updateBulkData();

    /**
     * @brief Updates reaction rates at steady state.
     * 
     * Computes ionization, attachment and effective reaction rates
     * both by counting the number of reaction occourred and 
     * by convolution.
     */
    void updateReactionRates();

    /**
     * @brief Assigns to each electron a collision type (or no collision).
     * 
     * By generating a random number for each electron and building the
     * collision matrix, a type of collision with a backgorund ga sparticle 
     * (or no collision) is assigned to every electron.
     */
    void updateCollisionMatrix();

    /**
     * @brief Performs electron collisions with the background gas.
     */
    void performCollisions();

    /**
     * @brief Checks if steady-state has been reached.
     */
    void checkSteadyState();

    /**
     * @brief Determines if the simulation should stop.
     * 
     * Determines if the simulation should stop because either
     * convergence criteria are fullfiled, or the maximum number of electrons
     * has been reached, or no electrons are left in the swarm,
     * or the maximum number of collisions have been reached.
     * 
     * @return True if stopping criteria are met.
     */
    bool endSimulation();

    /**
     * @brief Prints partial results on screen
     */
    void printOnScreen();

    /**
     * @brief Saves simulation results to file.
     * @param duration Duration of the simulation.
     * 
     * Saves in a file the simulation data and final results.
     * The duration of the simulation is expected to be counted in the main file. 
     */
    void saveResults(const int64_t duration) const;

    /**
     * @brief Returns the number of steady-state steps.
     * @return Number of steady-state steps.
     */
    const unsigned int get_count_sst() const { return count_sst; }

private:

    // Class Members:
    const unsigned  N0;     ///< initial electron population   
    const double N;         ///< gas density in m^-3
    
    const std::vector<std::string> gas; ///< subformulas of gas species
    std::vector<double> mgas;           ///< mass of gas species (in kg)
    std::vector<double> mix;            ///< fractions of individual species
            
    CrossSectionsData Xsec;         ///< cross-section data for gas species
                    
    const double w_err;             ///< error tolerance of mean energy (default: 1%)
    const double DN_err;            ///< error tolerance of diffusion constant (default: 1%)
    const unsigned Ne_max;          ///< maximum allowed number of electrons
    const unsigned long col_equ;    ///< minimum # of collisions for equilibrium (default: 20e6)
    const unsigned long col_max;    ///< max # of collisions (default: 20e6) 
    const bool conserve;            ///< force electron conservation (true), allow gain/loss (false)
    const bool isotropic;           ///< isotropic scattering (true) or not (false)

    const double W;                 ///< energy sharing in ionizing collision
    const double sqrtW;             ///< std::sqrt(W) -> pre computation
    const double sqrt1_W;           ///< std::sqrt(1-W) -> pre computation
    double E_max;                   ///< maximum electron energy
    const double EN;                ///< Electric field per-electron [Td]

    double T_sst = 0.0;             ///< equilibrium time
    unsigned int count_sst = 0;     ///< counts the number of steady state time steps
    int line = 1;                   ///< line number in output file:

    std::vector<double> t;          ///< vector of times [s]
    double dt;                      ///< current time step [s]
    double t_total = 0.0;           ///< sum of all times for all electrons [m]
    mc::MATRIX r;                   ///< 3D positions of electrons [m]
    mc::MATRIX v;                   ///< current velocity of electrons [m/s]
    std::array<double,3> a;         ///< acceleration of electrons [m/s^2] (const & uniform)
    mc::MATRIX v_int;               ///< time-integrated velocity
    mc::MATRIX v2_int;              ///< time-integrated squared velocity
    std::vector<double> v_abs;      ///< absolute value of velocities
    std::vector<double> E_in_eV;    ///< kinetic energies of electrons [eV]
    mc::MATRIX r_cations;           ///< 3D positions of cations; [m]
    mc::MATRIX r_anions;            ///< 3D positions of anions;  [m]

    CollisionData C;                ///< Collision indeces, Mass and Loss vectors
    unsigned long collisions = 0;   ///< total number collisions occourred

    std::vector<MeanData> mean;     ///< temporal mean data of electron swarm
    EnergyData E;                   ///< energy data
    BulkData bulk;                  ///< bulk transport data
    FluxData flux;                  ///< flux transport data
    RateDataConv rates_conv;        ///< reaction rates computed by convolution
    RateDataCount rates_count;      ///< reaction rates computed by counting

    unsigned int converge = 0;      ///< Status
            
    double EnergyLossElastic    = 0;    ///< Heat energy loss in elastic collisions
    double EnergyLossInelastic  = 0;    ///< Heat energy loss in inelastic collisions
    double EnergyLossIonization = 0;    ///< heat energy loss in ionization collisions

    std::default_random_engine gen;             ///< Random number generator
    std::uniform_real_distribution<> randu;     ///< Uniform distribution for random numbers
    std::normal_distribution<double> randn;     ///< Normal distribution for random numbers

    // Private Methods:

    /**
     * @brief Checks that the sum of gas fractions equals 1.
     * 
     * Called during initialization.
     * If the total fraction is not equal to 1, the fraction of
     * the last indexed component is adjusted accordingly. 
     */
    void checkFractionSum();

    /**
     * @brief Computes mass of the gas species in kg
     * 
     * Called during initialization.
     */
    void mass_in_kg();

    /**
     * @brief Sets initial position and velocity of electrons
     * 
     * Called during initialization
     */    
    void initialParticles(const std::array<double,3> & pos_xyz, const std::array<double,3> & sigma_xyz);

    /**
     * @brief Computes the cross product of two 3D arrays.
     * 
     * @param a First array.
     * @param b Second array.
     * @return Cross product a x b.
     */
    std::array<double, 3> cross_product(const std::array<double, 3>& a, const std::array<double, 3>& b) const;

    /**
     * @brief Performs elastic collisions for a given set of electrons.
     * 
     * @param ind Indices of the electrons.
     * @param Mass Masses of the gas-species undergoing elastic collsion.
     */
    void elasticCollision(const std::vector<size_t> & ind, const std::vector<double> & Mass);

    /**
     * @brief Performs inelastic collisions for a given set of electrons.
     * 
     * @param ind Indices of the electrons.
     * @param Loss Energy losses due to inelastic collisions.
     */
    void inelasticCollision(const std::vector<size_t> & ind, const std::vector<double> & Loss);

    /**
     * @brief Performs ionization collisions for a given set of electrons.
     * 
     * @param ind Indices of the electrons.
     * @param Loss Energy losses due to ionization.
     */
    void ionizationCollision(const std::vector<size_t> & ind, const std::vector<double> & Loss);

    /**
     * @brief Performs attachment collisions for a given set of electrons.
     * 
     * @param ind Indices of the electrons.
     */
    void attachmentCollision(const std::vector<size_t> & ind);
};

#endif  // MONTECARLO_HPP