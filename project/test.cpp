#include <iostream>
#include <vector>
#include <string>
#include <random>
#include "CrossSectionsData.hpp"
#include "CollisionData.hpp"
#include "Common.hpp"

// Dummy mc namespace and constants for compilation if not included elsewhere
// namespace mc { constexpr double q0 = 1.602e-19; constexpr double me = 9.109e-31; }

int main() {
    try {
        // Setup gas mixture and parameters
        std::vector<std::string> species = {"N2", "N", "O2", "O", "H2O", "H2"};
        std::vector<double> mix = {0.2, 0.2, 0.2, 0.2, 0.1, 0.1};
        double E_max = 1e5;
        double p = 101325; // Pressure in Pa
        double T = 300;    // Temperature in K

        double density = p / (mc::kB * T); // Density in m^-3

        // Construct CrossSectionsData
        CrossSectionsData xs(species, E_max, mix, density);

        // Print number and types of interactions for each specie
        for (size_t i = 0; i < species.size(); ++i) {
            std::vector<std::string> interaction_names;
            // Check all possible interaction types
            std::vector<std::string> types = {"EFFECTIVE", "IONIZATION", "ATTACHMENT", "EXCITATION", "ELASTIC"};
            for (const auto& type : types) {
                auto xs_tables = xs.get_Xsections(species[i], type);
                if (!xs_tables.empty()) {
                    interaction_names.push_back(type);
                }
            }
            std::cout << "Specie: " << species[i]
                      << " | Total interactions: " << interaction_names.size()
                      << " | Types: ";
            for (size_t j = 0; j < interaction_names.size(); ++j) {
                std::cout << interaction_names[j];
                if (j + 1 < interaction_names.size()) std::cout << ", ";
            }
            std::cout << std::endl;
        }


        // Prepare dummy mgas vector (masses for each species, just for test)
        std::vector<double> mgas(species.size(), 1.0);

        // Construct CollisionData
        CollisionData cd(xs, mgas);

        // Prepare dummy particle data
        int n_particles = 1e3;
        std::vector<double> E_in_eV(n_particles, 1e5); // All particles at 10 eV
        //std::vector<double> v_abs(n_particles, 1e6);    // All particles at 1e6 m/s

        // Test ComputeIndeces
        cd.ComputeIndeces(n_particles, xs, E_in_eV, mix, density);

        // Test getters
        std::cout << "Number of elastic collisions: " << cd.get_ind("ELASTIC").size() << std::endl;
        std::cout << "Number of excitation collisions: " << cd.get_ind("EXCITATION").size() << std::endl;
        std::cout << "Number of ionization collisions: " << cd.get_ind("IONIZATION").size() << std::endl;
        std::cout << "Number of attachment collisions: " << cd.get_ind("ATTACHMENT").size() << std::endl;
        std::cout << "Total collisions: " << cd.getCollisions() << std::endl;

        std::cout << "CollisionData test completed successfully." << std::endl;
    } catch (const std::exception& ex) {
        std::cerr << "Exception: " << ex.what() << std::endl;
        return 1;
    }
    return 0;
}