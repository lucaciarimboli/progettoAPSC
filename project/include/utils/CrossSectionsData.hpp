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

    // Constructor:
    CrossSectionsData(const std::vector<std::string> & species, const double E_max, 
                      const std::vector<double> & mix, const double N);

    // Getters:
    const std::vector<std::string>& get_gas() const { return gas; }
    const std::vector<double>& get_energy() const { return energy; }
    const double get_nu_max() const { return nu_max; }
    const std::vector<table>& get_full_xs_data() const { return Xsections; }
    const unsigned get_n_react() const;

    // Option 1: pass indexes of specie and interaction
    const std::vector<table> get_Xsections( const size_t specie, const mc::InteractionType interaction) const;
    // Option 2: pass strings with specie and interaction names
    const std::vector<table> get_Xsections( const std::string & specie, const std::string & interaction) const;

private:

    // Class Members:
    std::vector<std::string> gas;       // cell array of subformula of gas species
    std::vector<size_t> n_react;        // number of reactions for each specie
    std::vector<table> Xsections;       // Vector of cross-section data
    std::vector<double> energy;         // energy grid for cross sections
    double nu_max;                      // maximal collision frequency

    // Private methods:
    void count_react_and_fill_energy(const size_t i);
    void import_Xsec_data(const size_t offset, const size_t specie_indexe);
    void linear_interpolation(std::vector<double>& x, std::vector<double>& y, std::vector<double>& result);
    void maximalCollFreq(const std::vector<double> & mix, const double N);
};

#endif // CROSS_SECTIONS_DATA_HPP