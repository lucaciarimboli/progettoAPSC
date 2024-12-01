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
        for (size_t i = 1; i < energy_bins.size(); ++i) {
            if (energy_bins[i] <= energy_bins[i - 1]) {
            throw std::invalid_argument("Energy bins must be monotonically increasing.");
        }
        }
    };

    // Update energy statistics
    void update_energy(const std::vector<double>& E_in_eV, double dt, size_t ne) {

        E_sum += std::accumulate(E_in_eV.begin(), E_in_eV.end(), 0.0) * dt;

        // Count energy in every bin 
        for (const auto &en : E_in_eV) {
            // Loop over the bins except the last one
            for (size_t i = 0; i < energy.size()-1; i++) {
                if (en >= energy[i] && en < energy[i+1]) {
                    EEPF_sum[i]++;
                    break;
                }
            }
            // Last bin 
            if (en >= energy[energy.size() - 1]) {
                EEPF_sum.back()++;
            }
        }
    }

    // Normalize EEPF and compute EEDF
    void compute_distribution_function() {
        // Normalize EEPF
        double dE = energy[1] - energy[0];
        double total_sum = std::accumulate(EEPF_sum.begin(), EEPF_sum.end(), 0.0);
        for (size_t j = 0; j < energy.size(); ++j) {
            if(total_sum > 0.0){    //for safety
            EEPF[j] = EEPF_sum[j] / (total_sum * dE);
            }
        }
        
        // Compute EEDF
        for (size_t j = 0; j < energy.size(); ++j) {
            double energy_bin = energy[j];
            if (energy_bin == 0) {
                energy_bin = energy[1]; // Handle zero energy case
            }
            EEDF[j] = EEPF[j] / std::sqrt(energy_bin);
        }  
    }

    // Setters and getters for energy data
    void set_E_mean(const double & t_total){ E_mean = E_sum / t_total; }
    double get_E_mean() const { return E_mean; }
    const std::vector<double>& get_EEPF() const { return EEPF; }
    const std::vector<double>& get_EEDF() const { return EEDF; }

private:
    std::vector<double> energy;    // Energy bins
    std::vector<int> EEPF_sum;     // Histogram data
    std::vector<double> EEPF;      // Normalized Electron Energy Probability Function
    std::vector<double> EEDF;      // Electrons Energy Distribution Function
    double E_sum;                  // Sum of energies
    double E_mean;                 // Mean energy
};

#endif