#ifndef REACTION_RATES_HPP
#define REACTION_RATES_HPP

#include <vector>
#include <map>
#include <string>
#include <array>
#include <cmath>
#include "utils/MeanData.hpp"
#include "utils/EnergyData.hpp"
#include "utils/CrossSectionsData.hpp"

#include "Common.hpp"


class RateDataBase {
public:

    // Constructor:
    RateDataBase() { rates.fill(0.0); }

    // Virtual destructor:
    virtual ~RateDataBase() = default;

    // Virtual Method:
    virtual void computeRates() = 0;

    // Getter:
    double getRate(const mc::InteractionType key) const { return rates[key]; }

protected:

    // Protected member:
    std::array<double,mc::INTERACTIONS> rates; // Reaction rates

    // Protected method:
    const std::string inter_to_string(mc::InteractionType interaction) const;
};


// Calculates the reaction rates by regression of electron number vs time
class RateDataCount: public RateDataBase {
public:

    // THIS CLASS COMPUTES ONLY EFFECTIVE, IONIZATION AND ATTACHMENT RATES
    // RATES RELATED TO EXCITATION AND ELASTIC COLLISIONS ARE SIMPLY SET TO 0.0

    RateDataCount(const double & dens, bool cons) 
        : conserve(cons), N(dens)
    {
        rates_errors.fill(0.0);
    }

    // Public Method:
    void computeRates() override;

    // Setters:
    void setConserve(bool cons) { conserve = cons; }
    void setTime(const std::vector<double>& t,const unsigned int & count_sst);
    void setParticles(const std::vector<MeanData> & mean, const unsigned int & count_sst);

    // Getter:
    double get_errors(const mc::InteractionType rate_key) const { return rates_errors[rate_key]; }
    
private:

    // Private Members:
    std::array<double, mc::INTERACTIONS> rates_errors; // Errors in reaction rates

    bool conserve; // conserve (true) electron number after ionizatzion/attachment or not (false)
    double N;      // Gas number density in [m^-3]

    std::vector<double> x; // linear time interval
    std::array<std::vector<int>,mc::PARTICLES_TYPES> particles; // Number of particles of each type at each time of vector "x"

    // Private methods:
    void computeRate(const std::vector<double>& x, const std::vector<double>& y, const int & rate_key);
    void computeNonConserved();
    void computeConserved();

};


struct spec_rate {
    double rate;                      // Reaction rate
    std::vector<double> sigma;        // Cross-section data for the reaction
    size_t specie;                    // Specie (identified by its index in the mix)
    mc::InteractionType interaction;  // Type of interaction
    std::string reaction;             // Type of reaction
};


// Calculates the reaction rates by convolution of the electron number with a kernel
class RateDataConv : public RateDataBase {
    public:

    // Constructor:
    RateDataConv( const CrossSectionsData & xs, const EnergyData & en, const std::vector<double> & mix);

    // Public Method:
    void computeRates() override;

    // Getter:
    const std::vector<spec_rate> & getSpecificRates() const { return specific_rates; }

    private:

    // Class Members:
    std::vector<spec_rate> specific_rates; // Reaction Rates specific to each type, reaction and interaction.
    
    const CrossSectionsData& Xsec;                // cross_section data
    const EnergyData& E;                          // Energy data
    const std::vector<double>& mix;               // Mixture fractions

    // Private Methods:
    void linear_interpolation(const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& result) const;
    double convolution(const std::vector<double>& sigma, const std::vector<double>& EEPF, const std::vector<double>& sqrt_energy) const;
};

#endif // REACTION_RATES_HPP