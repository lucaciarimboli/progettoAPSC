#include "utils/EnergyData.hpp"

EnergyData::EnergyData(const double& E_max, const double& E_step):
    dE(E_step),
    E_sum(0.0),
    E_mean(0.0)
{
    const size_t n_bins = static_cast<size_t>(E_max / dE);
    if (n_bins < 2) {
        throw std::invalid_argument("Energy bins must have at least two elements.");
    }

    energy.reserve(n_bins);
    sqrt_E.reserve(n_bins);
    for (size_t i = 0; i < n_bins; i++){
        energy.push_back(i * dE);
        sqrt_E.push_back(std::sqrt(energy.back()));
    }
    // Note that energy.back() <= E_max

    EEPF_sum.resize(n_bins, 0);
    EEPF.resize(n_bins, 0.0);
    EEDF.resize(n_bins, 0.0);
};

void EnergyData::mean_energy(const mc::MATRIX& v2_int, const double& t_total){
    // Compute mean energy by time-averaging energies of all electron trajectories:

    std::vector<double> E_int(v2_int.size());
    std::transform(v2_int.begin(), v2_int.end(), E_int.begin(), [](const std::array<double, 3>& v2i){
        return 0.5 * mc::me * (v2i[0] + v2i[1] + v2i[2]) / mc::q0;
    });
    E_sum += std::accumulate(E_int.cbegin(), E_int.cend(), 0.0);
    E_mean = E_sum / t_total;
}

void EnergyData::compute_distribution_function(const std::vector<double>& E_in_eV) {

    // Count energy in every bin 
    for (auto it = E_in_eV.cbegin(); it != E_in_eV.cend(); it++) {
        const size_t bin_index = std::min(static_cast<size_t>(*it / dE), energy.size() - 1);
        EEPF_sum[bin_index]++;
    }

    // Normalize EEPF and compute EEDF
    const double total_sum = static_cast<double>(std::accumulate(EEPF_sum.cbegin(), EEPF_sum.cend(), 0));

    for (size_t j = 1; j < energy.size(); j++) {
        EEPF[j] = EEPF_sum[j] / (total_sum * dE);
        EEDF[j] = EEPF[j] / sqrt_E[j];
    }

    // Avoid division by 0 (energy[0] = 0.0):
    EEPF[0] = EEPF_sum[0] / (total_sum * dE);
    EEDF[0] = EEPF[0] / sqrt_E[1];
}