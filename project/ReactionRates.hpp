#ifndef REACTION_RATES_HPP
#define REACTION_RATES_HPP

#include <vector>
#include <map>
#include <string>
#include <array>
#include <cmath>
#include "MeanData.hpp"
#include "EnergyData.hpp"
#include "cross_s.h"

class cross_sect;
class MeanData;
class EnergyData;

class RateDataBase {
public:
    RateDataBase()
    {
        rates.fill(0.0);
    }

    virtual ~RateDataBase() = default;

    // Pure virtual function to compute rates:
    virtual void computeRates() = 0;

    // Getter for rates:
    double getRate(const std::string& key) const { return rates.at(key); }

protected:
    std::array<double,INTERACTIONS> rates; // Reaction rates
};

// Calculates the reaction rates by regression of electron number vs time
class RateDataCount: public RateDataBase {
public:

    // THIS CLASS COMPUTES ONLY EFFECTIVE, IONIZATION AND ATTACHMENT RATES
    // RATES RELATED TO EXCITATION AND ELASTIC COLLISIONS ARE SIMPLY SET TO 0.0

    RateDataCount(const double & dens, bool cons) 
        : N(dens), conserve(cons)
    {
        rates_errors.fill(0.0);
    }

    // Compute the reaction rates based on the time and particle data.
    void computeRates() override
    {
        if (!conserve) {
            computeNonConserved();
        } else {
            computeConserved();
        }
    }

    // Setters:
    void setConserve(bool cons) { conserve = cons; }
    void setTime(const std::vector<double>& t,const unsigned int & count_sst);
    void setParticles(const std::vector<MeanData> & mean, const unsigned int & count_sst);
    
private:

    std::array<double, INTERACTIONS> rates_errors; // Errors in reaction rates

    bool conserve; // conserve (1) electron number after ionizatzion/attachment or not (0)
    double N; // Gas number density in m^-3

    std::vector<double> x; // linear time interval
    std::array<std::vector<int>,PARTICLES_TYPES> particles; // Number of particles of each type at each time of vector "x"

    // Private methods:
    void computeRate(const std::vector<double>& x, const std::vector<double>& y, const int & rate_key);
    void computeNonConserved();
    void computeConserved();

};

struct spec_rate {
    double rate;              // Reaction rate
    int specie;               // Index of the specie
    int interaction;          // Type of interaction
    std::string reaction;     // Type of reaction
};

// Calculates the reaction rates by convolution of the electron number with a kernel
class RateDataConv : public RateDataBase {
    public:
    RateDataConv( const cross_sect & sigma, const std::vector<double> mixx, const EnergyData & en) 
        : Sigma(sigma), mix(mixx), E(en)
    {
        // Set the correct size for the specific rates vector:
        size_t specific_rates_size = 0;
        for( cross_sect & s : Sigma) {
            specific_rates_size += s.tab.size();
        }
        specific_rates.reserve(specific_rates_size);
    }

    void computeRates() override
    {
        // Iterator to specific rates vector:
        auto it = specific_rates.begin();

        // Loop over all species in the gas mixture:
        for(size_t specie = 0; specie < mix.length(); specie++) {

            // Iterate over all interactions:
            for( auto & table : Sigma(specie).tab) {

                // Compute the reaction rate for element "t"
                spec_rate rr;
                rr.rate = convolution(E, table.energy, table.sect);
                rr.specie = specie;
                rr.interaction = table.interact;
                rr.reaction = table.react;
                specific_rates.push_back(R);

                // Update total reaction rate
                rates[rr.interaction] += R.rate * mix[specie];
            }
        }
        // Compute the effective rate: 
        //rates[EFFECTIVE] = rates[IONIZATION] - rates[ATTACHMENT];
    }

    setSigma(const cross_sect & sigma) { Sigma = sigma; } // Set the cross-section data

    private:
    std::vector<spec_rate> specific_rates; // Reaction Rates specific to each type, reaction and interaction.
    std::vector<cross_sect> Sigma;         // cross_section data
    std::vector<double> mix;               // Fractions of individual species in the gas as a vector
    EnergyData E;                          // Energy data

    double convolution(const std::vector<double> & x, const std::vector<double> & y);
    std::vector<double> linear_interpolation(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& xq);
    // linear interpolation requires C++20 for std::lerp
};

#endif // REACTION_RATES_HPP