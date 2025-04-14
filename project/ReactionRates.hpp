#ifndef REACTION_RATES_HPP
#define REACTION_RATES_HPP

#include <vector>
#include <map>
#include <string>
#include <array>
#include "MeanData.hpp"


class RateDataBase {
public:
    RateDataBase()
    {
        rates["ion_tot"] = 0.0;
        rates["att_tot"] = 0.0;
        rates["eff"] = 0.0;
    }

    virtual ~RateDataBase() = default;

    // Pure virtual function to compute rates:
    virtual void computeRates() = 0;

    // Getter for rates:
    double getRate(const std::string& key) const { return rates.at(key); }

protected:
    std::map<std::string, double> rates; // Reaction rates
};

// Calculates the reaction rates by regression of electron number vs time
class RateDataCount: public RateDataBase {
public:

    RateDataCount(const double & dens, bool cons) 
        : N(dens), conserve(cons)
    {
        rates["ion_tot_err"] = 0.0;
        rates["att_tot_err"] = 0.0;
        rates["eff_err"] = 0.0;
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
    bool conserve; // conserve (1) electron number after ionizatzion/attachment or not (0)
    double N; // Gas number density in m^-3

    std::vector<double> x; // linear time interval
    std::array<std::vector<int>,PARTICLES_TYPES> particles; // Number of particles of each type at each time of vector "x"

    // Private methods:
    void computeRate(const std::vector<double>& x, const std::vector<double>& y, const std::string& rate_key, const std::string& err_key);
    void computeNonConserved();
    void computeConserved();

};

// Calculates the reaction rates by convolution of the electron number with a kernel
class RateDataConv : public RateDataBase {
    public:
    RateDataConv() 
        : RateDataBase(N)
    {
        rates["ela"] = 0.0;
        rates["ela_tot"] = 0.0;
        rates["ion"] = 0.0;
        rates["att"] = 0.0;
        rates["exc"] = 0.0;

    }

    void computeRates() override
    {
        // ...
    }     

    private:
    
};

#endif // REACTION_RATES_HPP