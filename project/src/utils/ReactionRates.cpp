#include "utils/ReactionRates.hpp"

void RateDataCount::setTime(const std::vector<double>& t, const unsigned int & count_sst){
    x.assign(t.cend() - count_sst, t.cend());
    std::transform(x.begin() + 1, x.end(), x.begin() + 1, [this](double xx) { return xx - x[0]; });
    x[0] = 0.0;
}

void RateDataCount::setParticles(const std::vector<MeanData> & mean, const unsigned int & count_sst) {
    //if (count_sst > mean.size()) {throw std::out_of_range("count_sst exceeds the size of mean");}

    for (size_t i = 0; i < mc::PARTICLES_TYPES; i++) {

        // Set the size for particle vector
        for (auto& p : particles) p.resize(count_sst);

        // Iterator to mean elements starting from sst
        auto it = mean.cend() - count_sst;

        for (size_t j = 0; j < count_sst; j++) {
            // Get the particle counts from MeanData:
            const std::array<int,mc::PARTICLES_TYPES>& mean_particles = (it + j)->get_particles();
    
            for (size_t k = 0; k < mc::PARTICLES_TYPES; k++)
                particles[k][j] = mean_particles[k]; // Fill the particles vector with the number of particles of each type
        }
    } 
}

void RateDataCount::computeRate(const std::vector<double>& x, const std::vector<double>& y, const int & rate_key) {
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
    
    rates[rate_key] = std::accumulate(ratio.begin(), ratio.end(), 0.0) / (ratio.size()) / N;
    rates_errors[rate_key] = std::sqrt(std::inner_product(ratio.begin(), ratio.end(), ratio.begin(), 0.0) / (ratio.size() - 1)) / (x.size() - 1) / N;
}

void RateDataCount::computeNonConserved() {

    // Compute effective ionization rate:
    std::vector<double> y(particles[mc::ELECTRONS].size());

    // Define y vector as logarithimcal electrons gain
    if(particles[mc::ELECTRONS][0] == 0) throw std::invalid_argument("Number of electrons cannot be zero");
    y[0] = 0.0;
    std::transform(particles[mc::ELECTRONS].begin() + 1, particles[mc::ELECTRONS].end(), y.begin() + 1, [this](int part_0i) {
        if( part_0i == 0) throw std::invalid_argument("Number of electrons cannot be zero");
        return std::log(static_cast<double>(part_0i) / static_cast<double>(particles[mc::ELECTRONS][0]));
    });

    computeRate(x, y, mc::EFFECTIVE);

    double nu_eff = rates[mc::EFFECTIVE] * N;
    std::transform(x.begin(), x.end(), x.begin(), [nu_eff, this](int xi) {
        return (std::exp(nu_eff * xi) - 1.0) / nu_eff * static_cast<double>(particles[mc::ELECTRONS][0]);;
    });

    // Compute ionization rate:
    std::transform(particles[mc::CATIONS].begin(), particles[mc::CATIONS].end(), y.begin(), [this](int part_ij) {
        return static_cast<double>(part_ij - particles[mc::CATIONS][0]);      // y elements here are int but they are casted to double for coherence with y definition
    });
    computeRate(x, y, mc::IONIZATION);

    // Compute attachment rate:
    std::transform(particles[mc::ANIONS].begin(), particles[mc::ANIONS].end(), y.begin(), [this](int part_ij) {
        return static_cast<double>(part_ij - particles[mc::ANIONS][0]);      // same as for ionization rates
    });
    computeRate(x, y, mc::ATTACHMENT);
}

void RateDataCount::computeConserved() {
    // y[0] is for effective ionization, y[1] for ionization and y[2] for attachment
    std::array<std::vector<double>,mc::PARTICLES_TYPES> y;
    size_t count_sst = particles[0].size();

    double initial_electrons = static_cast<double>(particles[mc::ELECTRONS][0]);
    if(initial_electrons == 0) throw std::invalid_argument("Number of electrons cannot be zero");

    // Define y vector as normalized particles gain
    for( size_t i = 0; i < mc::PARTICLES_TYPES; i++) {
        y[i].resize(count_sst);
        std::transform(particles[i].begin(), particles[i].end(), y[i].begin(), [this,i,initial_electrons](int part_ij) {
            return (static_cast<double>(part_ij - particles[i][0]) / initial_electrons);
        });
    }

    // Compute effective ionization, ionization and attachment rates:
    for (int i = 0; i < 3; i++){
        computeRate(x, y[i], i);
    }
    // Equivalent to:
    // computeRate(x, y[ELECTRONS], EFFECTIVE);
    // computeRate(x, y[CATIONS], IONIZATION);
    // computeRate(x, y[ANIONS], ATTACHMENT);
}

// Compute the convolution for vectors x,y over the energy grid provided in E
double RateDataConv::convolution(const std::vector<double>& x, const std::vector<double>& y) {

    // By construction in CrossSectionData.cpp,
    // the Xsections energy grid covers all the energy values of the simulation grid.
    // Convolution and linear interpolation can be called safely without further checks.

    // Get energy values and electrons energy probability function:
    const std::vector<double>& energy = E.get_energy();
    const std::vector<double>& EEPF = E.get_EEPF(); 
    double dx = energy[1] - energy[0]; // Assume uniform grid spacing
    double rate = 0.0;

    // Interpolate cross-section values to match the energy grid of E
    //const std::vector<double> sigma_f = linear_interpolation(x, y, energy);
    size_t i = 0; // Index for the x vector (energy grid for xs)
    for (size_t n = 0; n < energy.size(); n++) {

        // Interpolate xs value over the energy grid:
        const double xq = energy[n]; 
        const double t = (xq - x[i]) / (x[i + 1] - x[i]);
        const double sigma_f = y[i] + t * (y[i + 1] - y[i]);
        
        // Update index i:
        while( i < x.size() - 1 && xq > x[i + 1]) i++;

        // Update convolution integral result:
        rate += EEPF[n] * std::sqrt(xq) * sigma_f * dx;
    }

    // Return the computed rate:
    return std::sqrt(2.0 * mc::q0 / mc::me) * rate;   
}