#include <iostream>
#include <vector>
#include <string>
#include "CrossSectionsData.hpp"

double maximalCollFreq(const CrossSectionsData & Xsec, const std::vector<double> & mix);

int main() {
    try {
        // Define all allowed species
        std::vector<std::string> species = {"H2", "H2O", "N", "N2", "O", "O2"};

        // Define mock mixture fractions (must sum to 1.0)
        std::vector<double> mix = {0.2, 0.2, 0.2, 0.2, 0.1, 0.1};

        // Create a CrossSectionsData object
        CrossSectionsData crossSectionsData(species);

        // Print gas species
        std::cout << "Gas species: ";
        for (const auto& gas : crossSectionsData.get_gas()) {
            std::cout << gas << " ";
        }
        std::cout << std::endl;

        // Print energy grid size:
        std::cout << "Energy grid size: " << crossSectionsData.get_energy().size() << std::endl;

        // Print energy grid
        std::cout << "Energy grid: ";
        const auto& energy_grid = crossSectionsData.get_energy();
        size_t grid_size = energy_grid.size();
        for (size_t i = 0; i < std::min(grid_size, size_t(5)); ++i) {
            std::cout << energy_grid[i] << " ";
        }
        if (grid_size > 10) {
            std::cout << "... ";
        }
        for (size_t i = (grid_size > 5 ? grid_size - 5 : 5); i < grid_size; ++i) {
            std::cout << energy_grid[i] << " ";
        }
        std::cout << std::endl;

        // Print full cross-section data for each species
        const auto& full_xs_data = crossSectionsData.get_full_xs_data();
        for (size_t i = 0; i < species.size(); ++i) {
            std::cout << "Cross-section data for " << species[i] << ":" << std::endl;
            for (const auto& table : full_xs_data[i]) {
                std::cout << "  Reaction: " << table.react << std::endl;
                std::cout << "  Average energy: " << table.en_avg << " eV" << std::endl;
                std::cout << "  Interaction type: " << table.interact << std::endl;
                std::cout << "  Cross-section values: ";
                size_t section_size = table.section.size();
                for (size_t i = 0; i < std::min(section_size, size_t(3)); ++i) {
                    std::cout << table.section[i] << " ";
                }
                if (section_size > 6) {
                    std::cout << "... ";
                }
                for (size_t i = (section_size > 3 ? section_size - 3 : 3); i < section_size; ++i) {
                    std::cout << table.section[i] << " ";
                }
                std::cout << std::endl;
                std::cout << "  ------------------------" << std::endl;
            }
        }

        // Compute total cross-section
        std::vector<double> total_xsection = crossSectionsData.compute_total_Xsection(mix);

        // Print total cross-section
        std::cout << "\n\nTotal cross-section: ";
        size_t total_size = total_xsection.size();
        for (size_t i = 0; i < std::min(total_size, size_t(3)); ++i) {
            std::cout << total_xsection[i] << " ";
        }
        if (total_size > 6) {
            std::cout << "... ";
        }
        for (size_t i = (total_size > 3 ? total_size - 3 : 3); i < total_size; ++i) {
            std::cout << total_xsection[i] << " ";
        }
        std::cout << std::endl;

        // Calculate maximal collision frequency
        double nu_max = maximalCollFreq(crossSectionsData, mix);
        std::cout << "\n\nMaximal collision frequency: " << nu_max << " s^-1" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Exception occurred: " << e.what() << std::endl;
    }

    return 0;
}

double maximalCollFreq(const CrossSectionsData & Xsec, const std::vector<double> & mix) {
    // calculates maximal collision rate for gas mixture (in s^-1)
    double nu_max = 0.0;
    double N = 1.0; // Number density (example value, adjust as needed)
    double q0 = 1.6e-19; // Elementary charge in Coulombs
    double me = 9.11e-31; // Electron mass in kg

    // Extract total_cross section data from Xsec object:
    const std::vector<double> & energy = Xsec.get_energy();
    std::vector<double> sigma_tot = Xsec.compute_total_Xsection(mix);

    // Check if the size is coherent:
    if(energy.size() != sigma_tot.size()){
        std::cerr << "Error: energy and total_xsection vectors have different sizes!" << std::endl;
        return 0.0;
    }

    // Calculate the maximal collision frequency:
    for(size_t i = 0; i < energy.size(); i++){
        double nu = N * sigma_tot[i] * std::sqrt(2.0 * energy[i] * q0 / me); // s^-1
        if(nu > nu_max) nu_max = nu;
    }

    return nu_max;
}