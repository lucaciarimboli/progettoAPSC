#include "ReactionRates.hpp"

void RateDataCount::set_x(const std::vector<double>& t,const unsigned int & count_sst){
    x.assign(t.cend() - count_sst, t.cend());
    std::transform(x.begin() + 1, x.end(), x.begin() + 1, [this](double val) { return val - x[0]; });
    x[0] = 0.0;
}

void RateDataCount::set_particles(const std::vector<MeanData> & mean, const unsigned int & count_sst) {
    if (count_sst > mean.size()) {throw std::out_of_range("count_sst exceeds the size of mean");}
    for (size_t i = 0; i < PARTICLES_TYPES; ++i) {

        // Set the size for particle vector
        for (auto& yy : y) {
            yy.resize(count_sst);
        }

        // Iterator to mean elements starting from sst
        auto it = mean.cend() - count_sst;

        // Initial (at initial time x[0]) number of particles
        std::array<int, PARTICLES_TYPES> initial_particles = it->get_particles();
        initial_electron_population = initial_particles[ELECTRONS];

        for (size_t i = 0; i < count_sst; i++) {
            it++; // Move to the next element
            const auto & mean_particles = it->get_particles(); // Get the particle counts from MeanData
    
            for (size_t j = 0; j < PARTICLES_TYPES; ++j) {
                y[j][i] = mean_particles[j]; // Fill the particles vector with the number of particles of each type 
            }
        }
    } 
}

void computeRate(const std::vector<double>& x, const std::vector<double>& y, const std::string& rate_key, const std::string& err_key) {
    
    std::vector<double> ratio(x.size() - 1); // ratio is element-wise division of y by time x

    for (size_t i = 1; i < x.size(); ++i) {
        ratio[i - 1] = y[0][i] / x[i]; 
    }
    
    rates["eff"] = std::accumulate(ratio.begin(), ratio.end(), 0.0) / (ratio.size()) / N;
    rates["eff_err"] = std::sqrt(std::inner_product(ratio.begin(), ratio.end(), ratio.begin(), 0.0) / (ratio.size() - 1)) / (x.size() - 1) / N;
}

void RateDataCount::computeNonConserved() {

    // Compute effective ionization rate:

    // Define y vector as logarithimcal particles gain
    std::transform(y[0].begin(), y[0].end(), y[0].begin(), [initial_particles, j](double val) {
        return std::log(val) - std::log(initial_electron_population);
    });

    computeRate(x, y[0], "eff", "eff_err");

    double nu_eff = rates["eff"] * N;
    std::transform(x.begin(), x.end(), x.begin(), [nu_eff, initial_electron_population](double xi) {
        return (std::exp(nu_eff * xi) - 1.0) / nu_eff * initial_electron_population;
    });

    // Compute ionization and attachment rates:

    // Define y vector as particles gain:
    for( i = 1; i < PARTICLES_TYPES; i) {
        std::transform(y[i].begin(), y[i].end(), y[i].begin(), [initial_particles, j](double val) {
            return val - initial_particles[j];
        });
    }

    computeRate(x, y[1], "ion_tot", "ion_tot_err");
    computeRate(x, y[2], "att_tot", "att_tot_err");
}

void RateDataConv::computeConserved() {
    // Define y vector as normalized particles gain
    for( i = 1; i < PARTICLES_TYPES; i) {
        std::transform(y[i].begin(), y[i].end(), y[i].begin(), [initial_particles, j](double val) {
            return (val - initial_particles[j]) / inital_electron_population;
        });
    }

    // Compute effective ionization, ionization and attachment rates:
    computeRate(x, y[0], "eff", "eff_err");
    computeRate(x, y[1], "ion_tot", "ion_tot_err");
    computeRate(x, y[2], "att_tot", "att_tot_err");
}
