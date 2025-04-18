#include "ReactionRates.hpp"

void RateDataCount::setTime(const std::vector<double>& t, const unsigned int & count_sst){
    x.assign(t.cend() - count_sst, t.cend());
    std::transform(x.begin() + 1, x.end(), x.begin() + 1, [this](double xx) { return xx - x[0]; });
    x[0] = 0.0;
}

void RateDataCount::setParticles(const std::vector<MeanData> & mean, const unsigned int & count_sst) {
    if (count_sst > mean.size()) {throw std::out_of_range("count_sst exceeds the size of mean");}
    for (size_t i = 0; i < PARTICLES_TYPES; i++) {

        // Set the size for particle vector
        for (auto& p : particles) p.resize(count_sst);

        // Iterator to mean elements starting from sst
        auto it = mean.cend() - count_sst;

        for (size_t i = 0; i < count_sst; i++) {
            // Get the particle counts from MeanData:
            const std::array<int,PARTICLES_TYPES>& mean_particles = (it + i)->get_particles();
    
            for (size_t j = 0; j < PARTICLES_TYPES; j++)
                particles[j][i] = mean_particles[j]; // Fill the particles vector with the number of particles of each type
        }
    } 
}

void RateDataCount::computeRate(const std::vector<double>& x, const std::vector<double>& y, const std::string& rate_key, const std::string& err_key) {
    // Check correctness of x,y.
    if(x.size() < 2) {
        throw std::invalid_argument("x vector must have at least two elements");
    }
    if(y.size() != x.size()) {
        throw std::invalid_argument("x and y vectors must have the same size");
    }
    // Remark that inside MonteCarlo.cpp, if this method is called x,y have the same size >10.

    std::vector<double> ratio(x.size() - 1); // ratio is element-wise division of y by time x

    for (size_t i = 1; i < x.size(); i++) {
        ratio[i - 1] = y[i] / x[i]; 
    }
    
    rates["eff"] = std::accumulate(ratio.begin(), ratio.end(), 0.0) / (ratio.size()) / N;
    rates["eff_err"] = std::sqrt(std::inner_product(ratio.begin(), ratio.end(), ratio.begin(), 0.0) / (ratio.size() - 1)) / (x.size() - 1) / N;
}

void RateDataCount::computeNonConserved() {

    // Compute effective ionization rate:
    std::vector<double> y(particles[ELECTRONS].size());

    // Define y vector as logarithimcal electrons gain
    if(particles[ELECTRONS][0] == 0) throw std::invalid_argument("Number of electrons cannot be zero");
    y[0] = 0.0;
    std::transform(particles[ELECTRONS].begin() + 1, particles[ELECTRONS].end(), y.begin() + 1, [this](int part_0i) {
        if( part_0i == 0) throw std::invalid_argument("Number of electrons cannot be zero");
        return std::log(static_cast<double>(part_0i) / static_cast<double>(particles[ELECTRONS][0]));
    });

    computeRate(x, y, "eff", "eff_err");

    double nu_eff = rates["eff"] * N;
    std::transform(x.begin(), x.end(), x.begin(), [nu_eff, this](int xi) {
        return (std::exp(nu_eff * xi) - 1.0) / nu_eff * static_cast<double>(particles[ELECTRONS][0]);;
    });

    // Compute ionization rate:
    std::transform(particles[CATIONS].begin(), particles[CATIONS].end(), y.begin(), [this](int part_ij) {
        return static_cast<double>(part_ij - particles[CATIONS][0]);      // y elements here are int but they are casted to double for coherence with y definition
    });
    computeRate(x, y, "ion_tot", "ion_tot_err");

    // Compute attachment rate:
    std::transform(particles[ANIONS].begin(), particles[ANIONS].end(), y.begin(), [this](int part_ij) {
        return static_cast<double>(part_ij - particles[ANIONS][0]);      // same as for ionization rates
    });
    computeRate(x, y, "att_tot", "att_tot_err");
}

void RateDataCount::computeConserved() {
    // y[0] is for effective ionization, y[1] for ionization and y[2] for attachment
    std::array<std::vector<double>,PARTICLES_TYPES> y;
    size_t count_sst = particles[0].size();

    double initial_electrons = static_cast<double>(particles[ELECTRONS][0]);
    if(initial_electrons == 0) throw std::invalid_argument("Number of electrons cannot be zero");

    // Define y vector as normalized particles gain
    for( size_t i = 0; i < PARTICLES_TYPES; i++) {
        y[i].resize(count_sst);
        std::transform(particles[i].begin(), particles[i].end(), y[i].begin(), [this,i](int part_ij) {
            return (static_cast<double>(part_ij - particles[i][0]) / initial_electrons);
        });
    }

    // Compute effective ionization, ionization and attachment rates:
    computeRate(x, y[ELECTRONS], "eff", "eff_err");
    computeRate(x, y[CATIONS], "ion_tot", "ion_tot_err");
    computeRate(x, y[ANIONS], "att_tot", "att_tot_err");
}
