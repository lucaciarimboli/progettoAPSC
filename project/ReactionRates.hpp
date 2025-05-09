#ifndef REACTION_RATES_HPP
#define REACTION_RATES_HPP

#include <vector>
#include <map>
#include <string>
#include <array>
#include <cmath>
#include "MeanData.hpp"
#include "EnergyData.hpp"
#include "CrossSectionsData.hpp"

class CrossSectionsData;
class MeanData;
class EnergyData;

class RateDataBase {
public:
    RateDataBase() {
        rates.fill(0.0);
    }

    virtual ~RateDataBase() = default;

    // Pure virtual function to compute rates:
    virtual void computeRates() = 0;

    // Getter for rates:
    double getRate(const INTER& key) const { return rates[key]; }

protected:
    std::array<double,INTERACTIONS> rates; // Reaction rates

    // Converts INTER enum to string for printing:
    std::string INTER_to_string(const INTER & interaction) const {
        switch (interaction) {
            case EFFECTIVE: return "EFFECTIVE";
            case IONIZATION: return "IONIZATION";
            case ATTACHMENT: return "ATTACHMENT";
            case EXCITATION: return "EXCITATION";
            case ELASTIC: return "ELASTIC";
            default: return "UNKNOWN";
        }
    }
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

    // Getters:
    double get_errors(const INTER & rate_key) const { return rates_errors[rate_key]; }

    // Printer:
    void printRates() const
    {
        std::cout << "\nReaction Rates:" << std::endl;
        for (int i = 0; i < 3; i++) {
            std::cout << "Rate " << INTER_to_string(static_cast<INTER>(i)) << ": " << rates[i] << std::endl;
            std::cout << "Error: " << rates_errors[i] << std::endl;
            std::cout << "-----------" << std::endl;
        }
    }
    
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
    std::string specie;               // Index of the specie
    std::string interaction;          // Type of interaction
    std::string reaction;     // Type of reaction
};

// Calculates the reaction rates by convolution of the electron number with a kernel
class RateDataConv : public RateDataBase {
    public:
    RateDataConv( const CrossSectionsData & xs, const EnergyData & en) 
        : Xsec(xs), E(en)
    {
        // Set the correct size for the specific rates vector:
        size_t specific_rates_size = 0;
        const std::vector<std::vector<table>> & Sigma = Xsec.get_full_xs_data();
        for( const std::vector<table> & t : Sigma) {
            specific_rates_size += t.size();
        }
        specific_rates.reserve(specific_rates_size);
    }

    void computeRates() override
    {
        // Iterator to specific rates vector:
        auto it = specific_rates.begin();

        // Vector containing mixture fractions:
        std::vector<double> mix = Xsec.get_mix();
        // Vector containing the energy values:
        std::vector<double> energy_grid = Xsec.get_energy();

        // Loop over all species in the gas mixture:
        for(size_t specie = 0; specie < mix.size(); specie++) {

            const std::vector<table> & Sigma = Xsec.get_full_xs_data()[specie];

            // Iterate over all interactions:
            for( const table & t : Sigma) {

                // Compute the reaction rate for element "t"
                spec_rate rr;
                rr.rate = convolution(energy_grid, t.section);
                rr.specie = Xsec.get_gas()[specie];
                rr.interaction = INTER_to_string(t.interact);
                rr.reaction = t.react;
                specific_rates.push_back(rr);

                // Update total reaction rate
                rates[t.interact] += rr.rate * mix[specie];
            }
        }
    }

    // Setters:
    void setSigma(const CrossSectionsData & xs) { Xsec = xs; } // Set the cross-section data
    void setEnergy(const EnergyData & en) { E = en; } // Set the energy data
    // Getters:
    const std::vector<spec_rate> & getSpecificRates() const { return specific_rates; }   // Get the specific rates
    // Printers:
    void printRates() const
    {
        std::cout << "\nReaction Rates:" << std::endl;
        for (int i = 0; i < 5; i++) {
            std::cout << "Rate " << INTER_to_string(static_cast<INTER>(i)) << ": " << rates[i] << std::endl;
            std::cout << "-----------" << std::endl;
        }
    }
    
    void printSpecificRates() const {
        std::cout << "\nSpecific Reaction Rates:" << std::endl;
        for (const auto& rr : specific_rates) {
            std::cout << "Specie: " << rr.specie << std::endl;
            std::cout << "Interaction: " << rr.interaction << std::endl;
            std::cout << "Reaction: " << rr.reaction << std::endl;
            std::cout << "Rate: " << rr.rate << std::endl;
            std::cout << "-----------" << std::endl;
        }
    }

    private:
    std::vector<spec_rate> specific_rates; // Reaction Rates specific to each type, reaction and interaction.
    CrossSectionsData Xsec;                // cross_section data
    EnergyData E;                          // Energy data

    double convolution(std::vector<double> x, std::vector<double> y);
    std::vector<double> linear_interpolation(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& xq);
};

#endif // REACTION_RATES_HPP