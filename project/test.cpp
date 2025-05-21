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
        std::cout << "Mass vector: ";
        for (auto m : cd.getMass()) std::cout << m << " ";
        std::cout << "\nLoss vector: ";
        for (auto l : cd.getLoss()) std::cout << l << " ";
        std::cout << "\nElastic collision indices: ";
        for (auto i : cd.get_ind("ELASTIC")) std::cout << i << " ";
        std::cout << "\nExcitation collision indices: ";
        for (auto i : cd.get_ind("EXCITATION")) std::cout << i << " ";
        std::cout << "\nIonization collision indices: ";
        for (auto i : cd.get_ind("IONIZATION")) std::cout << i << " ";
        std::cout << "\nAttachment collision indices: ";
        for (auto i : cd.get_ind("ATTACHMENT")) std::cout << i << " ";
        std::cout << "\nTotal collisions: " << cd.getCollisions() << std::endl;

        std::cout << "CollisionData test completed successfully." << std::endl;
    } catch (const std::exception& ex) {
        std::cerr << "Exception: " << ex.what() << std::endl;
        return 1;
    }
    return 0;
}