#ifndef ENERGYDATA_H
#define ENERGYDATA_H

#include <vector>
#include <cmath>
#include <numeric>

class EnergyData {
public:
    // Default constructor
    EnergyData() = default;
    // Initialization constructor
    EnergyData(const std::vector<double>& energy_bins)
        : energy(energy_bins),
          EEPF_sum(energy_bins.size(), 0), 
          EEPF(energy_bins.size(), 0.0),
          EEDF(energy_bins.size(), 0.0),
          E_sum(0.0),
          E_mean(0.0)
    {
        // Check that "energy" format is correct:
        if (energy_bins.size() < 2) {
            throw std::invalid_argument("Energy bins must have at least two elements.");
        }
        /*for (size_t i = 1; i < energy_bins.size(); ++i) {
            if (energy_bins[i] <= energy_bins[i - 1]) {
            throw std::invalid_argument("Energy bins must be monotonically increasing.");
        }
        }*/
    };

    // Update energy statistics
    void update_energy(const std::vector<double>& E_in_eV, const double dt, const unsigned int ne, const double t_total) {

        E_sum += std::accumulate(E_in_eV.begin(), E_in_eV.end(), 0.0) * dt;
        E_mean = E_sum / t_total;

        // Count energy in every bin 
        for (unsigned int i = 0; i < ne; i++) {
            const double en = E_in_eV[i];
            // Use binary search to find the appropriate bin
            if (en >= energy.front() && en <= energy.back()) {
                const double step_size = energy[1] - energy[0];
                const size_t bin_index = std::min(static_cast<size_t>(en / step_size), energy.size() - 1);
                EEPF_sum[bin_index]++;
            }
        }
    }

    // Normalize EEPF and compute EEDF
    void compute_distribution_function() {
        // Normalize EEPF
        double dE = energy[1] - energy[0];
        double total_sum = std::accumulate(EEPF_sum.begin(), EEPF_sum.end(), 0.0);
        if(total_sum > 0.0){

            for (size_t j = 1; j < energy.size(); j++) {
                EEPF[j] = EEPF_sum[j] / (total_sum * dE);
                EEDF[j] = EEPF[j] / std::sqrt(energy[j]);
            }

            // Handle separately the energy level E=0 to avoid division by 0:
            EEPF[0] = EEPF_sum[0] / (total_sum * dE);
            EEDF[0] = EEPF[0] / std::sqrt(energy[1]);
        }
    }

    // Setters:
    // void set_energy(const std::vector<double>& energy_bins){ energy = energy_bins; }
    // Getters:
    const std::vector<double>& get_energy() const { return energy; }
    const double get_E_mean() const { return E_mean; }
    const std::vector<double>& get_EEPF() const { return EEPF; }
    const std::vector<double>& get_EEDF() const { return EEDF; }

private:
    const std::vector<double> energy;    // Energy bins
    std::vector<int> EEPF_sum;     // Histogram data
    std::vector<double> EEPF;      // Normalized Electron Energy Probability Function
    std::vector<double> EEDF;      // Electrons Energy Distribution Function
    double E_sum;                  // Sum of energies
    double E_mean;                 // Mean energy
};

#endif