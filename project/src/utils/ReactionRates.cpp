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
    if (!conserve) {
        computeNonConserved();
    } else {
        computeConserved();
    }
}

void RateDataCount::setTime(const std::vector<double>& t, const unsigned int & count_sst){
    x.assign(t.cend() - count_sst, t.cend());
    std::transform(x.begin() + 1, x.end(), x.begin() + 1, [this](double xx) { return xx - x[0]; });
    x[0] = 0.0;
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

        // Get the particle counts from MeanData:
        const std::array<int,mc::PARTICLES_TYPES>& mean_particles = it->get_particles();
    
        // Fill the particles vector with the number of particles of each type
        particles[mc::ELECTRONS].push_back(mean_particles[mc::ELECTRONS]);
        particles[mc::ANIONS].push_back(mean_particles[mc::ANIONS]);
        particles[mc::CATIONS].push_back(mean_particles[mc::CATIONS]);
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
    size_t count_sst = particles[mc::ELECTRONS].size();
    y[0].resize(count_sst, 0.0);
    y[1].resize(count_sst, 0.0);
    y[2].resize(count_sst, 0.0);

    const double initial_electrons = static_cast<double>(particles[mc::ELECTRONS][0]);
    //if(initial_electrons == 0) throw std::invalid_argument("Number of electrons cannot be zero");

    // Compute effective ionization, ionization and attachment rates:
    for (int i = 0; i < 3; i++){

        // Define y vector as normalized particles gain
        std::transform(particles[i].cbegin(), particles[i].cend(), y[i].begin(), [this,i,initial_electrons](int part_ij) {
            return (static_cast<double>(part_ij - particles[i][0]) / initial_electrons);
        });

        // Compute effective ionization, ionization and attachment rates:
        computeRate(x, y[i], i);
    }
}

RateDataConv::RateDataConv( const CrossSectionsData & xs, const EnergyData & en, const std::vector<double> & mix) 
    : Xsec(xs), E(en), mix(mix)
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
    for(spec_rate & rr : specific_rates) {

        // Compute the reaction rate for element "rr"
        rr.rate = convolution(rr.sigma);

        // Update total reaction rate
        rates[rr.interaction] += rr.rate * mix[rr.specie];
    }

    rates[mc::EFFECTIVE] = rates[mc::IONIZATION] - rates[mc::ATTACHMENT];
}

double RateDataConv::convolution(const std::vector<double>& sigma) const
{

    // Compute the rate by convolution of the cross section
    // with the probability density function of having an electron at a given energy level.

    // Get energy values and electrons energy probability function:
    const std::vector<double>& energy = E.get_energy();
    const std::vector<double>& EEPF = E.get_EEPF(); 
    const double& dx = E.get_dE(); // Assume uniform grid spacing
    double rate = 0.0;

    for (size_t j = 0; j < energy.size(); j++) {
        rate += EEPF[j] * std::sqrt(energy[j]) * sigma[j] * dx;
    }

    // Return the computed rate:
    return std::sqrt(2.0 * mc::q0 / mc::me) * rate;   
}