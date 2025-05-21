#ifndef CROSS_SECTIONS_Data_HPP
#define CROSS_SECTIONS_Data_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <filesystem>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <numeric>

#include "Common.hpp"

struct table {
    std::vector<double> section;
    size_t specie_index;
    double en_avg;
    std::string react;
    mc::InteractionType interact;
};

class CrossSectionsData {
public:

    // Constructor
    CrossSectionsData(const std::vector<std::string> & species, const double E_max, 
                      const std::vector<double> & mix, const double N);

    // Getters:
    const std::vector<std::string>& get_gas() const { return gas; }
    const std::vector<double>& get_energy() const { return energy; }
    const double get_nu_max() const { return nu_max; }
    const std::vector<table>& get_full_xs_data() const { return Xsections; }
    const unsigned get_n_react() const {
        unsigned sum = 0;
        for (auto n : n_react) sum += n;
        return sum - gas.size(); // Exclude effective xs data for each specie
    }
    // getters for specific cross-section data per specie-interaction type:
    // Option 1: pass indexes of specie and interaction
    const std::vector<table> get_Xsections( const size_t specie, const mc::InteractionType interaction) const;
    // Option 2: pass strings with specie and interaction names
    const std::vector<table> get_Xsections( const std::string & specie, const std::string & interaction) const;

private:
    // Members:
    std::vector<std::string> gas;       // cell array of subformula of gas species
    std::vector<size_t> n_react;        // number of reactions for each specie
    std::vector<table> Xsections;       // Vector of cross-section data
    std::vector<double> energy;         // energy grid for cross sections
    double nu_max;                      // maximal collision frequency

    // Private methods:

    // Counts the number of reactions for each specie and the number of cross-sections data points
    void count_react_and_fill_energy(const size_t i);

    // Imports cross-section data from .txt files
    void import_Xsec_data(const size_t offset, const std::string & current_specie);

    // Linear interpolation function
    void linear_interpolation(std::vector<double>& x, std::vector<double>& y, std::vector<double>& result);

    // Compute the maximal collision frequency:
    void maximalCollFreq(const std::vector<double> & mix, const double N);
};

#endif // CROSS_SECTIONS_DATA_HPP