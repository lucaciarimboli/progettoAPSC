#include <iostream>
#include <vector>
#include "CollisionData.hpp"

int main() {
    try {
        // Number of electrons and reactions
        const int n_particles = 5;
        const int n_react = 3;

        // Electron energies and velocities (mock data)
        std::vector<double> E_in_eV = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::vector<double> v_abs = {1e5, 1.5e5, 2e5, 2.5e5, 3e5};

        // Cross-section energy grid and cross-section data (mock data)
        std::vector<double> XS_energy = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::vector<std::vector<table>> XS(n_react, std::vector<table>(1));
        XS[0][0] = { {1e-20, 2e-20, 3e-20, 4e-20, 5e-20}, 0.0, "Elastic", ELASTIC };
        XS[1][0] = { {1e-21, 2e-21, 3e-21, 4e-21, 5e-21}, 10.0, "Excitation", EXCITATION };
        XS[2][0] = { {1e-22, 2e-22, 3e-22, 4e-22, 5e-22}, 15.0, "Ionization", IONIZATION };

        // Mixture fractions and masses (mock data)
        std::vector<double> mix = {0.5, 0.3, 0.2};
        std::vector<double> mgas = {2.0, 18.0, 28.0};

        // Other parameters
        double E_max = 5.0;
        double nu_max = 1e6;
        double density = 1e25;

        // Construct CollisionData object
        CollisionData collisionData(n_particles, n_react);

        // Build the collision matrix and compute indices
        collisionData.collisionmatrix(n_particles, E_in_eV, v_abs, XS_energy, XS, mix, mgas, E_max, nu_max, density);

        // Print collision matrix
        std::cout << "Collision Matrix (C):" << std::endl;
        for (int i = 0; i < n_react; ++i) {
            for (int j = 0; j < n_particles; ++j) {
                std::cout << collisionData.getC(i, j) << " ";
            }
            std::cout << std::endl;
        }

        // Print Mass and Loss vectors
        std::cout << "Mass vector: ";
        for (auto m : collisionData.getMass()) std::cout << m << " ";
        std::cout << std::endl;

        std::cout << "Loss vector: ";
        for (auto l : collisionData.getLoss()) std::cout << l << " ";
        std::cout << std::endl;

        // Print collision indices
        std::cout << "Elastic collision indices: ";
        for (auto idx : collisionData.get_ind("ELASTIC")) std::cout << idx << " ";
        std::cout << std::endl;

        std::cout << "Excitation collision indices: ";
        for (auto idx : collisionData.get_ind("EXCITATION")) std::cout << idx << " ";
        std::cout << std::endl;

        std::cout << "Ionization collision indices: ";
        for (auto idx : collisionData.get_ind("IONIZATION")) std::cout << idx << " ";
        std::cout << std::endl;

        std::cout << "Attachment collision indices: ";
        for (auto idx : collisionData.get_ind("ATTACHMENT")) std::cout << idx << " ";
        std::cout << std::endl;

        std::cout << "Total number of collisions: " << collisionData.getCollisions() << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Exception occurred: " << e.what() << std::endl;
    }
    return 0;
}