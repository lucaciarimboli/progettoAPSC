#include <iostream>
#include "ReactionRates.hpp"
#include "MeanData.hpp"

int main() {
    try {
        // Define constants
        const double gas_density = 1.0e25; // Example gas number density in m^-3
        const bool conserve = true;      // Example: do not conserve electron number
        const unsigned int count_sst = 10; // Example steady-state count

        // Create a RateDataCount object
        RateDataCount rateData(gas_density, conserve);

        // Define time vector
        std::vector<double> time = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0};

        // Set the time data
        rateData.setTime(time, count_sst);

        // Define MeanData objects
        std::vector<MeanData> meanData;
        for (size_t i = 0; i < time.size(); i++) {
            std::array<double, 3> position = {i * 1.0, i * 2.0, i * 3.0};
            std::array<double, 3> sigma = {0.1, 0.2, 0.3};
            std::array<int, PARTICLES_TYPES> particles = {
                100 - static_cast<int>(i * 3), // Electrons
                50 + static_cast<int>(i * 2),  // Cations
                30 + static_cast<int>(i * 1)   // Anions
            };
            meanData.emplace_back(particles, MATRIX{}, MATRIX{}, particles[ELECTRONS]);
        }

        // Set the particle data
        rateData.setParticles(meanData, count_sst);

        // Compute the reaction rates
        rateData.computeRates();

        // Print the computed rates
        std::cout << "Computed Rates:" << std::endl;
        std::cout << "Ionization Rate: " << rateData.getRate("ion_tot") << std::endl;
        std::cout << "Attachment Rate: " << rateData.getRate("att_tot") << std::endl;
        std::cout << "Effective Rate: " << rateData.getRate("eff") << std::endl;

        // Print the computed errors
        std::cout << "Ionization Rate Error: " << rateData.getRate("ion_tot_err") << std::endl;
        std::cout << "Attachment Rate Error: " << rateData.getRate("att_tot_err") << std::endl;
        std::cout << "Effective Rate Error: " << rateData.getRate("eff_err") << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Exception occurred: " << e.what() << std::endl;
    }

    return 0;
}