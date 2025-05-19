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

#include "Common.hpp"

struct table {
    std::vector<double> section;
    double en_avg;
    std::string react;
    mc::InteractionType interact;
};

struct table_tool {   // Just for temporary storage in the constructor
    std::vector<double> energy_tool;
    std::vector<double> section_tool;
};

class CrossSectionsData {
public:

    // Constructor
    CrossSectionsData(const std::vector<std::string> & species): gas(species)
    {
        // Check if the gas species are the same as the cross-section data
        std::array<std::string, 6> valid_species = {"H2", "H2O", "N", "N2", "O", "O2"};
        for(const auto & specie : gas) {
            if (std::find(std::begin(valid_species), std::end(valid_species), specie) == std::end(valid_species)) {
                throw std::invalid_argument("Invalid gas species: " + specie);
            }
        }

        // Set Xsections size:
        Xsections.resize(gas.size()); // Initialize Xsections with the correct size
        // VA BENE QUESTA OPERAZIONE O NO? CONSIDERANDO CHE Xsections E' UN VETTORE DI VETTORI?

        n_react = 0; // Initialize number of reactions

        // Xsec_tool is needed just to store temporarily the energy and cross-section values
        // for each specie and for each interaction before interpolating over a common energy grid
        std::vector<std::vector<table_tool>> Xsec_tool(gas.size());

        for( size_t i = 0; i < gas.size(); i++) {
            // import Xsec data from .txt files
            import_Xsec_data(i,Xsec_tool);
        }

        // Sort and deduplicate the energy vector to "extract" a grid shared among all the xsections
        sort_deduplicate(energy);
        // Ensure that the energy vector contains the energy level E = 0.0:
        if (energy[0] > 0.0) {
                 energy.insert(energy.begin(), 0.0);
        }
        // Adapt cross-section values to the common energy grid (by interpolation)
        fill_Xsections(Xsec_tool);
    }

    // Getters:
    const std::vector<std::string>& get_gas() const { return gas; }
    const std::vector<double>& get_energy() const { return energy; }
    const std::vector<std::vector<table>>& get_full_xs_data() const { return Xsections; }
    const int get_n_react() const { return n_react; }
    // getters for specific cross-section data per specie-interaction type:
    // Option 1: pass indexes of specie and interaction
    std::vector<table> get_Xsections( const size_t specie, const mc::InteractionType interaction) const {
        // Check if the specie and interaction indices are valid
        if (specie >= gas.size() || interaction >= Xsections[specie].size()) {
            throw std::out_of_range("Invalid specie or interaction index");
        }
        // Return the cross-section data for the specified specie and interaction
        std::vector<table> result;
        for (size_t i = 0; i < Xsections[specie].size(); i++) {
            if (Xsections[specie][i].interact == interaction) {
            result.push_back(Xsections[specie][i]);
            }
        }
        return result;
    }
    // Option 2: pass strings with specie and interaction names
    std::vector<table> get_Xsections( const std::string & specie, const std::string & interaction) const {
        size_t species_index = std::distance(gas.cbegin(), std::find(gas.cbegin(), gas.cend(), specie));

        // Check if the specie exists
        if (species_index == gas.size()) {
            throw std::invalid_argument("Invalid specie: " + specie);
        }

        // Map interaction strings to their corresponding enum values
        const std::map<std::string, mc::InteractionType> interaction_map = {
            {"EFFECTIVE", mc::EFFECTIVE},
            {"IONIZATION", mc::IONIZATION},
            {"ATTACHMENT", mc::ATTACHMENT},
            {"EXCITATION", mc::EXCITATION},
            {"ELASTIC", mc::ELASTIC}
        };

        // Check if the interaction exists
        auto it = interaction_map.find(interaction);
        if (it == interaction_map.end()) {
            throw std::invalid_argument("Invalid interaction: " + interaction);
        }
        // Get the cross-section data for the specified specie and interaction
        return get_Xsections(species_index, it->second);
    }
    
    // Computes the total cross-section
    std::vector<double> compute_total_Xsection(const std::vector<double> & mix) const {
        // Check if the mix vector is the same size as the gas vector
        if (mix.size() != gas.size()) {
            throw std::invalid_argument("Mix vector size must match the number of gas species");
        }

        std::vector<double> total_xsection(energy.size(), 0.0);

        // Iterate over all gas species
        for (size_t i = 0; i < gas.size(); i++) {
            // Iterate over all interactions
            for (size_t j = 0; j < Xsections[i].size(); j++) {

                // Multiply the cross-section by the corresponding mix fraction
                std::transform(total_xsection.begin(), total_xsection.end(),
                               Xsections[i][j].section.begin(),
                               total_xsection.begin(),
                               [mix, i](double tot_sigma, double sigma) {
                                   return tot_sigma + mix[i] * sigma;
                               });
            }
        }
        
        return total_xsection;
    }

private:
    // Members:
    std::vector<std::string> gas;              // cell array of subformula of gas species

    std::vector<std::vector<table>> Xsections; // Vector of cross-section data (each element is a specie)
    std::vector<double> energy;                // energy grid for cross sections
    int n_react;                               // total amount of possible reactions

    // Private methods:

    // Imports cross-section data from .txt files
    void import_Xsec_data(const size_t i, std::vector<std::vector<table_tool>> & Xsec_tool){

        // Path to the cross-section data file
        std::string path = "./Xsec/" + gas[i] + "/" + gas[i] + ".txt";

        // Map interactions indexes with their actual name in .txt data files
        const std::map<std::string, mc::InteractionType> int_map = {
            {"EFFECTIVE", mc::EFFECTIVE},
            {"IONIZATION", mc::IONIZATION},
            {"ATTACHMENT", mc::ATTACHMENT},
            {"EXCITATION", mc::EXCITATION},
            {"ELASTIC", mc::ELASTIC}
        };

        std::ifstream file;               
        std::string line;                  
        std::string avg_energy,formula;
        bool fl=0;
        unsigned n=0;
        unsigned counter=0;

        // Count the number of interactions in the file
        file.open(path,std::ifstream::in);
        if(file.is_open()){
            while(getline(file,line))
            {
                if (line.compare("IONIZATION")==0 ||
                    line.compare("ATTACHMENT")==0 || line.compare("EXCITATION")==0 ||
                    line.compare("ELASTIC") == 0){
                        n_react++;  // counts overall reactions (effective excluded)
                        n++;        // counts reactions of the current specie to reserve Xsections[i] memory
                    }
                
                else if(line.compare("EFFECTIVE")==0) n++;
            }
        }
        file.close();

        // Initialize i-esim vector of Xsections
        Xsections[i].resize(n);
        Xsec_tool[i].resize(n); 

        // Import data:
        file.open(path,std::ifstream::in);
        if(file.is_open()){
            int j=0;
            while(getline(file,line)) {
                // Set interaction type and average energy:
                if (line.compare("EFFECTIVE")==0 || line.compare("IONIZATION")==0 ||
                    line.compare("ATTACHMENT")==0 || line.compare("EXCITATION")==0 ||
                    line.compare("ELASTIC") == 0){

                    fl=1;
                    auto it = int_map.find(line);

                    // Set interaction type
                    Xsections[i][j].interact = it->second;
                
                    // Set reaction formula:
                    getline(file,formula);
                    Xsections[i][j].react=formula;
                
                    // Set average energy:
                    getline(file,avg_energy);
                    std::stringstream ss_en(avg_energy);
                    ss_en >> Xsections[i][j].en_avg;
                }

                // Set energy and cross sectional values:
                if(fl==1) {
                    if(line[0]=='-') counter++;
                    std::stringstream ss(line);
                    double x,y;  // x = energy level; y = cross section value

                    if(counter < 2) {
                        if(ss >> x >> y) {
                            // Store energy and cross-section values in the temporary vectors
                            Xsec_tool[i][j].energy_tool.push_back(x);
                            Xsec_tool[i][j].section_tool.push_back(y);

                            // Store energy values in the main energy vector
                            energy.push_back(x);
                        }
                    } else {
                        fl=0;
                        counter=0;
                        j++;
                    }
                }
            }
        }
    }

    // Sorts the input vector and removes duplicate elements
    void sort_deduplicate(std::vector<double>& v) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end(),
            [](double a, double b){
                return std::fabs(a - b) < 1e-10;  // tolerate small differences
            }), v.end());
    }

    // Fills missing energy values by interpolation to adapt the cross-sections to the "merged" energy grid
    void fill_Xsections(const std::vector<std::vector<table_tool>> & Xsec_tool){
        for(size_t i = 0; i < Xsec_tool.size(); i++) {
            for(size_t j = 0; j < Xsec_tool[i].size(); j++) {

                std::vector<double> e = Xsec_tool[i][j].energy_tool;
                std::vector<double> sigma = Xsec_tool[i][j].section_tool;
                double E_max = energy.back();

                // Add energy levels E=0 and E=E_max to fit the data for interpolation
                if( e[0] > 0.0) {
                    e.insert(e.begin(), 0.0);
                    sigma.insert(sigma.begin(), sigma.front());
                }
                if( e.back() < E_max) {
                    e.push_back(E_max);
                    sigma.push_back(sigma.back());
                }
                
                // Interpolate cross-section 
                Xsections[i][j].section = linear_interpolation(e, sigma, energy);
            }
        }
    }

    // Linear interpolation function
    std::vector<double> linear_interpolation(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& xq){

        if (x.empty() || y.empty() || x.size() != y.size()) {
           throw std::invalid_argument("Input vectors x and y must be non-empty and of the same size.");
        }

        std::vector<double> result(xq.size());

        for (size_t n = 0; n < xq.size(); n++) {
            double query = xq[n];
            size_t i = 0;

            if (query <= x.front()) {
            i = 0;
            } else if (query >= x.back()) {
            i = x.size() - 2;
            } else {
            while (i < x.size() - 1 && query > x[i + 1]) i++;
            }

            double t = (query - x[i]) / (x[i + 1] - x[i]);
            result[n] = y[i] + t * (y[i + 1] - y[i]);
        }

        return result;
    } 
};

#endif // CROSS_SECTIONS_DATA_HPP