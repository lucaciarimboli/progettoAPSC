#ifndef CROSS_SECTIONS_Data_HPP
#define CROSS_SECTIONS_Data_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <filesystem>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
#include <numeric>

//#include "Common.hpp"
#include "../Common.hpp" // for testing


/**
 * Struct that holds a single reaction cross section data.
 */
struct table {
    std::vector<double> section;
    size_t specie_index;
    double en_avg;
    std::string react;
    mc::InteractionType interact;
};



/**
 * @brief Struct that holds a single cross-section table.
 * 
 * Cross section data for the specified species are imported from .txt data files
 * that must be in the LXCat format.
 * A shared energy grid is defined after all the data have been imported and the cross
 * sections are interpolated on such grid.  
 */
class CrossSectionsData {
public:

    /**
     * @brief Constructor.
     * @param species List of gas species (subcomponents of the gas mixture).
     * @param E_max Maximum energy value for the energy grid.
     * @param mix Gas species mass fractions in the mix.
     * @param N Gas number density.
     */
    CrossSectionsData(const std::vector<std::string> & species, const double& E_max, 
                      const std::vector<double> & mix, const double& N);

    // Getters:

    /** 
     * @brief Returns list of gas species. 
     */
    const std::vector<std::string>& get_gas() const { return gas; }

    /**
     *  @brief Returns the energy grid.
     */
    const std::vector<double>& get_energy() const { return energy; }

    /**
     *  @brief Returns the maximum collision frequency.
     */
    const double get_nu_max() const { return nu_max; }

    /**
     *  @brief Returns the full cross-section dataset.
     */
    const std::vector<table>& get_full_xs_data() const { return Xsections; }

    /**
     *  @brief Returns total number of reactions across all species.
     */
    const unsigned get_n_react() const;

    /**
     * @brief Returns cross sections for a given species and interaction type.
     * @param specie Index of the gas species.
     * @param interaction Interaction type.
     */
    const std::vector<table> get_Xsections( const size_t& specie, const mc::InteractionType& interaction) const;

    /**
     * @brief Returns cross sections for a given species and interaction type.
     * @param specie Name of the gas species.
     * @param interaction Name of the interaction type.
     */
    const std::vector<table> get_Xsections( const std::string & specie, const std::string & interaction) const;
    
    /**
     *  @brief Removes effective cross sections from dataset.
     */
    void remove_effective_xs();

private:

    // Class Members:
    std::vector<std::string> gas;       ///< cell array of subformula of gas species
    std::vector<size_t> n_react;        ///< number of reactions for each specie
    std::vector<table> Xsections;       ///< Vector of cross-section data
    std::vector<double> energy;         ///< energy grid for cross sections
    double nu_max;                      ///< maximal collision frequency

    // Private methods:

    /**
     * @brief Used in constructor to allocate memory and unify the energy grid.
     * @param i Index of the species.
     * 
     * Counts the number of reactions for each specie and the number of cross-sections data points
     * and fills the energy vector with the energy levels from each reaction, sorted and deduplicated.
     */
    void count_react_and_fill_energy(const size_t& i);

    /**
     * @brief Imports cross-section data from LXCat file.
     * @param offset Index offset for current species.
     * @param specie_index Index of the current species.
     */
    void import_Xsec_data(const size_t& offset, const size_t& specie_index);

    /**
     * @brief Performs linear interpolation.
     * @param x Input x values (energy).
     * @param y Input y values (cross section values).
     * @param result Output interpolated values on energy grid.
     */
    void linear_interpolation(std::vector<double>& x, std::vector<double>& y, std::vector<double>& result);

    /**
     * @brief Computes maximal collision frequency.
     * @param mix Mixture fractions.
     * @param N Number density.
     */
    void maximalCollFreq(const std::vector<double> & mix, const double& N);
};

#endif // CROSS_SECTIONS_DATA_HPP