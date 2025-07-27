#include "utils/ReactionRates.hpp"

const std::string RateDataBase::inter_to_string(mc::InteractionType interaction) const {
    // Converts IntractionType to string
    switch (interaction) {
        case mc::EFFECTIVE:  return "EFFECTIVE";
        case mc::IONIZATION: return "IONIZATION";
        case mc::ATTACHMENT: return "ATTACHMENT";
        case mc::EXCITATION: return "EXCITATION";
        case mc::ELASTIC:    return "ELASTIC";
        default:             return "UNKNOWN";
    }
}

// Compute the reaction rates based on the time and particle data.
void RateDataCount::computeRates()
{
    if (conserve) {
        computeConserved();
    } else {
        computeNonConserved();
    }
}

void RateDataCount::setTime(const std::vector<double>& t, const unsigned int & count_sst){
    x.assign(t.cend() - count_sst, t.cend());
    const double x0 = x[0];
    std::transform(x.begin(), x.end(), x.begin(), [x0](double& xx) { return xx - x0; });
}

void RateDataCount::setParticles(const std::vector<MeanData> & mean, const unsigned int & count_sst) {
    //if (count_sst > mean.size()) {throw std::out_of_range("count_sst exceeds the size of mean");}

    // Set the size for particle vector
    particles[mc::ELECTRONS].clear();
    particles[mc::ANIONS].clear();
    particles[mc::CATIONS].clear();

    particles[mc::ELECTRONS].reserve(count_sst);
    particles[mc::ANIONS].reserve(count_sst);
    particles[mc::CATIONS].reserve(count_sst);

    for (auto it = mean.cend() - count_sst; it != mean.cend(); it++) {

        // Get the number of particles from MeanData:
        const std::array<int,mc::PARTICLES_TYPES>& p = it->get_particles();
    
        // Fill the particles vector with the number of particles of each type
        particles[mc::ELECTRONS].push_back(p[mc::ELECTRONS]);
        particles[mc::ANIONS].push_back(p[mc::ANIONS]);
        particles[mc::CATIONS].push_back(p[mc::CATIONS]);
    } 
}

void RateDataCount::computeRate(const std::vector<double>& x, const std::vector<double>& y, const int & rate_key) {

    std::vector<double> ratio; // ratio is element-wise division of y by time x
    const size_t ratio_sz = x.size() - 1;
    ratio.reserve(ratio_sz);

    for (size_t i = 1; i < x.size(); i++) {
        ratio.push_back(y[i] / x[i]);
    }
    
    // Compute mean and standard error of y(2:end)./x(2:end)
    const double mean = std::accumulate(ratio.begin(), ratio.end(), 0.0) / ratio_sz;
    const double std = std::sqrt(std::accumulate(ratio.begin(), ratio.end(), 0.0,
        [mean](double sum, double val) { return sum + (val - mean) * (val - mean);}
        ) / (ratio_sz - 1));
    rates[rate_key] = mean / N;
    rates_errors[rate_key] = std / std::sqrt(static_cast<double>(ratio_sz)) / N;
}

void RateDataCount::computeNonConserved() {

    // Define y vector as logarithimcal electrons gain:
    std::vector<double> y(particles[mc::ELECTRONS].size(), 0.0);
    const double log_ne0 = std::log(static_cast<double>(particles[mc::ELECTRONS][0]));
    std::transform(particles[mc::ELECTRONS].begin() + 1, particles[mc::ELECTRONS].end(), y.begin() + 1, [log_ne0](int part_0i) {
        // if( part_0i == 0) throw std::invalid_argument("Number of electrons cannot be zero");
        return std::log(static_cast<double>(part_0i)) - log_ne0;
    });

    // Compute effective rate:
    computeRate(x, y, mc::EFFECTIVE);
    const double ne0 = static_cast<double>(particles[mc::ELECTRONS][0]);
    const double nu_eff = rates[mc::EFFECTIVE] * N;
    std::transform(x.begin(), x.end(), x.begin(), [nu_eff, ne0](int xi) {
        return (std::exp(nu_eff * xi) - 1.0) / nu_eff * ne0;
    });

    // Compute ionization rate:
    const int nc0 = particles[mc::CATIONS][0];
    std::transform(particles[mc::CATIONS].begin(), particles[mc::CATIONS].end(), y.begin(), [nc0](int part_ij) {
        return static_cast<double>(part_ij - nc0);
    });
    computeRate(x, y, mc::IONIZATION);

    // Compute attachment rate:
    const int na0 = particles[mc::ANIONS][0];
    std::transform(particles[mc::ANIONS].begin(), particles[mc::ANIONS].end(), y.begin(), [na0](int part_ij) {
        return static_cast<double>(part_ij - na0);
    });
    computeRate(x, y, mc::ATTACHMENT);
}

void RateDataCount::computeConserved() {

    const double initial_electrons = static_cast<double>(particles[mc::ELECTRONS][0]);
    const size_t count_sst = particles[mc::ELECTRONS].size();

    std::vector<double> y(count_sst,0.0);

    // Compute effective ionization, ionization and attachment rates:
    for (int i = 0; i < 3; i++){

        // Define y vector as normalized particles gain
        const int part_i0 = particles[i][0];
        std::transform(particles[i].cbegin(), particles[i].cend(), y.begin(), [part_i0,initial_electrons](const int& part_ij) {
            return (static_cast<double>(part_ij - part_i0) / initial_electrons);
        });

        // Compute effective ionization, ionization and attachment rates:
        computeRate(x, y, i);
    }
}

RateDataConv::RateDataConv( const CrossSectionsData & xs, const EnergyData & en, const std::vector<double> & mix, const double& dE) 
    : Xsec(xs), E(en), mix(mix), factor(std::sqrt(2.0 * mc::q0 / mc::me) * dE)
{
    // Set the correct size for the specific rates vector:
    // ( "+ mix.size()" to account for effective xs data )
    specific_rates.reserve(Xsec.get_n_react() + mix.size());

    // Vector containing the energy values:
    const std::vector<double> & energy_grid_XS = Xsec.get_energy();

    // Loop over all reactions:
    for(const table & t : Xsec.get_full_xs_data()) {
        // Compute the reaction rate for element "t"
        spec_rate rr;
        rr.sigma.reserve(E.get_energy().size());
        linear_interpolation(energy_grid_XS, t.section, rr.sigma);
        rr.specie = t.specie_index;
        rr.interaction = t.interact;
        rr.reaction = t.react;
        specific_rates.push_back(rr);
    }    
}

void RateDataConv::linear_interpolation(const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& result) const
{

    // x is an energy grid and y is the corresponding xs data relative to one element.
    // The aim is to interpolate the values of y over the grid "energy".

    // By construction in CrossSectionData.cpp,
    // the Xsections energy grid covers all the energy values of the simulation grid.
    // Convolution and linear interpolation can be called safely without further checks.

    // x goes from 0.0 to E_max.
    // energy goes from 0.0 to the least multiple of dE lower than E_max
    // --> "x" range is equal or larger than "E_max" range

    // Get energy values and electrons energy probability function:
    const std::vector<double>& energy = E.get_energy();

    // Interpolate cross-section values to match the energy grid of E
    //const std::vector<double> sigma_f = linear_interpolation(x, y, energy);
    size_t i = 0; // Index for the x vector (energy grid for xs)
    for (auto it = energy.cbegin(); it != energy.cend(); it++) {

        const double xq = *it; // current energy grid node

        // Update index i:
        while( i < x.size() - 1 && xq > x[i + 1]) i++;

        // Interpolate xs value over the energy grid: 
        const double t = (xq - x[i]) / (x[i + 1] - x[i]);
        result.push_back(y[i] + t * (y[i + 1] - y[i]));
    }
}

void RateDataConv::computeRates(){

    const std::vector<double>& EEPF = E.get_EEPF(); 
    const std::vector<double>& sqrt_energy = E.get_sqrt_E();

    for(spec_rate & rr : specific_rates) {

        // Compute the reaction rate for element "rr":
        const std::vector<double>& sig = rr.sigma;
        double convolution = 0.0;
        for (size_t j = 0; j < sqrt_energy.size(); j++) {
            convolution += EEPF[j] * sqrt_energy[j] * sig[j];
        } // convolution integral
        //rr.rate = factor * convolution(rr.sigma, EEPF, sqrt_energy);
        rr.rate = factor * convolution;

        // Update total reaction rate
        rates[rr.interaction] += rr.rate * mix[rr.specie];
    }

    rates[mc::EFFECTIVE] = rates[mc::IONIZATION] - rates[mc::ATTACHMENT];
}

/*
double RateDataConv::convolution(const std::vector<double>& sigma, const std::vector<double>& EEPF,
      const std::vector<double>& sqrt_energy) const
{
    // Compute the rate by convolution of the cross section
    // with the probability density function of having an electron at a given energy level.
    
    double rate = 0.0;
    for (size_t j = 0; j < sqrt_energy.size(); j++) {
        rate += EEPF[j] * sqrt_energy[j] * sigma[j];
    }

    return rate;
}
*/