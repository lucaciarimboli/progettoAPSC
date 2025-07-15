#ifndef ENERGYDATA_H
#define ENERGYDATA_H

#include <vector>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <algorithm>

#include "Common.hpp"

class EnergyData {
public:
    // Constructors:
    EnergyData() = default;
    EnergyData(const double& E_max, const double& E_step);

    // Public methods:
    void mean_energy(const mc::MATRIX& v2_int, const double& t_total);
    void energy_bins(const std::vector<double>& E_in_eV);
    void compute_distribution_function();

    // Getters:
    const std::vector<double>& get_energy() const { return energy; }
    const double& get_dE() const { return dE; }
    const double& get_E_mean() const { return E_mean; }
    const std::vector<double>& get_EEPF() const { return EEPF; }
    const std::vector<double>& get_EEDF() const { return EEDF; }

private:

    // Class members:
    std::vector<double> energy;    // Energy bins
    const double dE;                     // Energy step
    std::vector<int> EEPF_sum;           // Histogram data
    std::vector<double> EEPF;            // Normalized Electron Energy Probability Function
    std::vector<double> EEDF;            // Electrons Energy Distribution Function
    double E_sum;                        // Sum of energies
    double E_mean;                       // Mean energy
};

#endif