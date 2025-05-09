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
    std::vector<double> y(particles[ELECTRONS].size());

    // Define y vector as logarithimcal electrons gain
    if(particles[ELECTRONS][0] == 0) throw std::invalid_argument("Number of electrons cannot be zero");
    y[0] = 0.0;
    std::transform(particles[ELECTRONS].begin() + 1, particles[ELECTRONS].end(), y.begin() + 1, [this](int part_0i) {
        if( part_0i == 0) throw std::invalid_argument("Number of electrons cannot be zero");
        return std::log(static_cast<double>(part_0i) / static_cast<double>(particles[ELECTRONS][0]));
    });

    computeRate(x, y, EFFECTIVE);

    double nu_eff = rates[EFFECTIVE] * N;
    std::transform(x.begin(), x.end(), x.begin(), [nu_eff, this](int xi) {
        return (std::exp(nu_eff * xi) - 1.0) / nu_eff * static_cast<double>(particles[ELECTRONS][0]);;
    });

    // Compute ionization rate:
    std::transform(particles[CATIONS].begin(), particles[CATIONS].end(), y.begin(), [this](int part_ij) {
        return static_cast<double>(part_ij - particles[CATIONS][0]);      // y elements here are int but they are casted to double for coherence with y definition
    });
    computeRate(x, y, IONIZATION);

    // Compute attachment rate:
    std::transform(particles[ANIONS].begin(), particles[ANIONS].end(), y.begin(), [this](int part_ij) {
        return static_cast<double>(part_ij - particles[ANIONS][0]);      // same as for ionization rates
    });
    computeRate(x, y, ATTACHMENT);
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
double RateDataConv::convolution(std::vector<double> x, std::vector<double> y) {
    // Check correctness of x,y
    if(y.size() != x.size()) {
        throw std::invalid_argument("x and y vectors must have the same size");
    }

    // Ensure energy starts at 0
    if (x[0] > 0) {
        y.insert(y.begin(), y[0]);
        x.insert(x.begin(), 0.0);
    }

    // Add a large energy value (1e10 eV) at the end
    y.push_back(y.back());
    x.push_back(1e10);

    // Get energy values and electrons energy probability function:
    const std::vector<double>& energy = E.get_energy();
    const std::vector<double>& EEPF = E.get_EEPF(); 

    // Interpolate cross-section values to match the energy grid of E
    std::vector<double> sigma_f = linear_interpolation(x, y, energy);

    double dx = energy[1] - energy[0]; // Assume uniform grid spacing
    double rate = 0.0;

    // Perform the convolution integral:
    for (size_t i = 0; i < EEPF.size(); i++) {
        rate += EEPF[i] * std::sqrt(energy[i]) * sigma_f[i] * dx;
    }

    // Return the computed rate:
    double me = 9.10938291e-31; // electron mass
    double q0 = 1.60217657e-19; // electron charge
    return std::sqrt(2.0 * q0 / me) * rate;   
}

// Performs linear interpolation
std::vector<double> RateDataConv::linear_interpolation(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& xq){

    if (x.empty() || y.empty() || x.size() != y.size()) {
        throw std::invalid_argument("Input vectors x and y must be non-empty and of the same size.");
    }

    std::vector<double> result(xq.size());

    for (size_t n = 0; n < xq.size(); n++) {
        double query = xq[n];
        size_t i = 0;

        if (query <= x.front()) {
        i = 0;
        } else if (query >= x.back()) {
        i = x.size() - 2;
        } else {
        while (i < x.size() - 1 && query > x[i + 1]) i++;
        }

        double t = (query - x[i]) / (x[i + 1] - x[i]);
        result[n] = y[i] + t * (y[i + 1] - y[i]);
    }

    return result;
} 