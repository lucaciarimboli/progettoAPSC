#ifndef ENERGYDATA_H
#define ENERGYDATA_H

#include <vector>
#include <cmath>
#include <numeric>
#include <stdexcept>

class EnergyData {
public:
    // Constructors:
    EnergyData() = default;
    EnergyData(const std::vector<double>& energy_bins);

    // Public methods:
    void update_energy(const std::vector<double>& E_in_eV, const double dt, const unsigned int ne, const double t_total);
    void compute_distribution_function();

    // Getters:
    const std::vector<double>& get_energy() const { return energy; }
    const double get_E_mean() const { return E_mean; }
    const std::vector<double>& get_EEPF() const { return EEPF; }
    const std::vector<double>& get_EEDF() const { return EEDF; }

private:

    // Class members:
    const std::vector<double> energy;    // Energy bins
    std::vector<int> EEPF_sum;           // Histogram data
    std::vector<double> EEPF;            // Normalized Electron Energy Probability Function
    std::vector<double> EEDF;            // Electrons Energy Distribution Function
    double E_sum;                        // Sum of energies
    double E_mean;                       // Mean energy
};

#endif