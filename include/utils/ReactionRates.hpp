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

#include <execution> 

/**
 * @class RateDataBase
 * @brief Abstract base class for reaction rates.
 *
 * This class defines the common interface and shared data structure
 * for the two different reaction rates computation strategies (convolution and counting).
 */
class RateDataBase {
public:

    /**
     *  @brief Constructor.
     * 
     */
    RateDataBase() { rates.fill(0.0); }

    /**
     *  @brief Virtual destructor.
     * 
     */
    virtual ~RateDataBase() = default;

    /**
     * @brief Computes the reaction rates.
     *
     * Pure virtual method to be implemented by derived classes.
     */
    virtual void computeRates() = 0;

    /**
     * @brief Returns the computed rate for a specific interaction.
     * @param key Type of interaction.
     * @return Reaction rate corresponding to the given key.
     */
    double getRate(const mc::InteractionType key) const { return rates[key]; }

protected:

    // Protected member:
    std::array<double,mc::INTERACTIONS> rates; ///< Reaction rates

    // Protected method:

    /**
     * @brief Converts an interaction enum to string for labeling/debugging.
     * @param interaction The interaction type.
     * @return A human-readable string representation.
     */
    const std::string inter_to_string(mc::InteractionType interaction) const;
};


/**
 * @class RateDataCount
 * @brief Computes reaction rates from particle count regression over time.
 *
 * This class computes effective, ionization and attachment rates by
 * conting the number of occourrences.
 */
class RateDataCount: public RateDataBase {
public:

    /**
     * @brief Constructor.
     * @param dens Gas number density
     * @param cons Flag for fixed electron number.
     * @param N0 Initial number of electrons.
     */
    RateDataCount(const double & dens, bool cons, const unsigned& N0) 
        : conserve(cons), N(dens), initial_electrons(static_cast<double>(N0))
    {
        rates_errors.fill(0.0);
    }

    // Public Method:
    /**
     *  @brief Computes reaction rates by counting.
     */
    void computeRates() override;

    // Setters:

    /**
     * @brief Sets conservation flag.
     * @param cons True if electrons number is conserved, false otherwise.
     */
    void setConserve(bool & cons) { conserve = cons; }
    /**
     * @brief Sets the time vector (starting at steady state).
     * @param time Time steps.
     * @param count_sst Number of steady state time steps.
     */
    void setTime(const std::vector<double>& time,const size_t & count_sst);
    void setParticles(const std::vector<MeanData> & mean, const size_t & count_sst);

    // Getter:
    /**
     * @brief Returns the error associated with a reaction rate.
     * @param rate_key Type of interaction.
     * @return Relative error of the computed rate.
     */
    double get_errors(const mc::InteractionType rate_key) const { return rates_errors[rate_key]; }
    
private:

    // Private Members:
    std::array<double, mc::INTERACTIONS> rates_errors; ///< Errors in reaction rates

    bool conserve; ///< conserve (true) electron number after ionizatzion/attachment or not (false)
    double N;      ///< Gas number density in [m^-3]

    std::vector<double> t;       ///< time vector starting from 0 after steady state
    std::vector<int> electrons;  ///< number of electrons at each time step
    std::vector<int> anions;     ///< number of anions at each time step
    std::vector<int> cations;    ///< number of cations at each time step
    double initial_electrons;    ///< Initial e population casted to double

    // Private methods:

    /**
     * @brief Compute the mean reaction rate and its standard error from time-series data.
     * 
     * @param x Vector of times (starting at steady state)
     * @param y Vector of normalized particles gain values (e.g. electrons gain, cations gain...).
     * @param rate_key The key indicating which reaction rate is being computed.
     */
    void computeRate(const std::vector<double>& x, const std::vector<double>& y, const mc::InteractionType& rate_key);
    /**
     *  @brief Computes rate when not enforcing electron conservation.
     */
    void computeNonConserved();
    /**
     *  @brief Computes rate when enforcing electron conservation.
     */
    void computeConserved();
};

/**
 * @struct spec_rate
 * @brief Structure containing details of a specific reaction rate.
 * 
 * This structure is needed in RateDataConv to store the reaction rates
 * specific to every species in the mix.
 */
struct spec_rate {
    double rate;                      // Reaction rate
    std::vector<double> sigma;        // Cross-section data for the reaction
    size_t specie;                    // Specie (identified by its index in the mix)
    mc::InteractionType interaction;  // Type of interaction
    std::string reaction;             // Type of reaction
};


/**
 * @class RateDataConv
 * @brief Computes reaction rates by energy convolution using cross-section data.
 *
 * Uses a convolution of the energy distribution and cross sections
 * to compute reaction rates for different species and interactions.
 * It also stores the reaction rates relative to every single specie in the mix.
 */
class RateDataConv : public RateDataBase {
    public:

    /**
     * @brief Constructor.
     * @param xs Cross section database.
     * @param en Electrons energy distribution data.
     * @param mix Mixture fractions.
     * @param dE Energy step.
     */
    RateDataConv( const CrossSectionsData & xs, const EnergyData & en, const std::vector<double> & mix, const double& dE);

    /**
     *  @brief Computes all reaction rates via convolution.
     */
    void computeRates() override;

    /**
     * @brief Returns the list of specific reaction rates.
     * @return Vector of spec_rate structures.
     */
    const std::vector<spec_rate> & getSpecificRates() const { return specific_rates; }

    private:

    // Class Members:
    std::vector<spec_rate> specific_rates; // Reaction Rates specific to each type, reaction and interaction.
    
    const CrossSectionsData& Xsec;                ///< cross_section data
    const EnergyData& E;                          ///< Energy data
    const std::vector<double>& mix;               ///< Mixture fractions

    const double factor;                          ///< Pre-computed "sqrt(2*q0/me)*dE" for performance

    // Private Methods:
    /**
     * @brief Performs linear interpolation of cross-section data.
     * @param x X-values (energies).
     * @param y Y-values (cross-sections).
     * @param result Interpolated result vector.
 */
    void linear_interpolation(const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& result) const;
};

#endif // REACTION_RATES_HPP