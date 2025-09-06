#ifndef ENERGYDATA_HPP
#define ENERGYDATA_HPP

#include <vector>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <algorithm>

#include "Common.hpp"

/**
 * @class EnergyData
 * @brief Computes and stores electron energy distribution data.
 *
 * This class processes the kinetic energies of electrons obtained 
 * from a Monte Carlo simulation. 
 * It contains the energy grid for the simulation and it provides methods to
 * compute the mean energy and to construct both the Electron Energy Distribution 
 * Function (EEDF) and the Electron Energy Probability Function (EEPF).
 */
class EnergyData {
public:
    
    /**
     * @brief Default constructor.
     */
    EnergyData() = default;
    /**
     * @brief Constructor.
     * @param E_max Maximum energy considered [eV].
     * @param E_step Energy step width [eV].
     * 
     * Builds the energy grid over which computations are performed in the simulation. 
     */
    EnergyData(const double& E_max, const double& E_step);

    // Public methods:

    /**
     * @brief Computes the time-averaged mean energy of electrons.
     * 
     * @param v2_int Integral over time of velocity squared for all electrons.
     * @param t_total Total time of simulation for all electrons
     */
    void mean_energy(const mc::MATRIX& v2_int, const double& t_total);

    /**
     * @brief Computes the energy distribution functions.
     * @param E_in_eV Electron energies [eV].
     *
     * Builds a histogram of the electron energy distribution and normalizes it 
     * to obtain the EEPF and EEDF.
     */
    void compute_distribution_function(const std::vector<double>& E_in_eV);

    // Getters:
    // Getters:
    /**
     *  @brief Get the energy grid.
     *  @return Energy grid.
     */
    const std::vector<double>& get_energy() const { return energy; }

    /**
     *  @brief Get the element-wise square root of the energy grid.
     *  @return element-wise square root of energy.
     */
    const std::vector<double>& get_sqrt_E() const { return sqrt_E; }

    /**
     *  @brief Get the energy step.
     *  @return Energy step.
     */
    const double& get_dE() const { return dE; }

    /**
     *  @brief Get the time-averaged mean energy.
     *  @return Time-averaged mean energy.
     */
    const double& get_E_mean() const { return E_mean; }

    /**
     *  @brief Get the electron energy probability function.
     *  @return Electron energy probability function.
     */
    const std::vector<double>& get_EEPF() const { return EEPF; }

    /**
     *  @brief Get the electron energy distribution function.
     *  @return Electron energy distribution function.
     */
    const std::vector<double>& get_EEDF() const { return EEDF; }

private:

    // Class members:
    std::vector<double> energy;          ///< Energy grid
    std::vector<double> sqrt_E;          ///< Square root of energy values
    const double dE;                     ///< Energy step
    std::vector<int> EEPF_sum;           ///< Histogram data
    std::vector<double> EEPF;            ///< Normalized Electron Energy Probability Function
    std::vector<double> EEDF;            ///< Electrons Energy Distribution Function
    double E_sum;                        ///< Sum of energies
    double E_mean;                       ///< Mean energy
};

#endif  // ENERGYDATA_HPP