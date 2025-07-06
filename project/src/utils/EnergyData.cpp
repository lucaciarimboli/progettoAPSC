#include "utils/EnergyData.hpp"

EnergyData::EnergyData(const std::vector<double>& energy_bins):
    energy(energy_bins),
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

void EnergyData::mean_energy(const mc::MATRIX& v2_int, const double& t_total){
    // Compute mean energy by time-averaging energies of all electronm trajectories:

    std::vector<double> E_int(v2_int.size());
    std::transform(v2_int.begin(), v2_int.end(), E_int.begin(), [](const std::array<double, 3>& v2i){
        return 0.5 * mc::me * (v2i[0] + v2i[1] + v2i[2]) / mc::q0;
    });
    E_sum += std::accumulate(E_int.cbegin(), E_int.cend(), 0.0);
    E_mean = E_sum / t_total;
}

void EnergyData::energy_bins(const std::vector<double>& E_in_eV) {
    // Count energy in every bin 
    for (unsigned int i = 0; i < E_in_eV.size(); i++) {
        const double en = E_in_eV[i];
        // Use binary search to find the appropriate bin
        if (en >= energy.front() && en <= energy.back()) {
            const double step_size = energy[1] - energy[0];
            const size_t bin_index = std::min(static_cast<size_t>(en / step_size), energy.size() - 1);
            EEPF_sum[bin_index]++;
        }
    }
}

void EnergyData::compute_distribution_function() {
    // Normalize EEPF and compute EEDF
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