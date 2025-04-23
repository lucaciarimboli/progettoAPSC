#include <iostream>
#include <vector>
#include "Interpolation.hpp"

int main() {
    try {
        // Define input data points (x and y)
        std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
        std::vector<double> y = {0.0, 2.0, 4.0, 6.0, 8.0};

        // Define query points (xq) where we want to interpolate
        std::vector<double> xq = {0.5, 1.5, 2.5, 3.5};

        // Create an Interpolation object
        Interpolation interpolator;

        // Perform linear interpolation
        std::vector<double> yq = interpolator.interpolation(x, y, xq);

        // Print the results
        std::cout << "Linear Interpolation Results:" << std::endl;
        std::cout << "xq\t\tyq" << std::endl;
        for (size_t i = 0; i < xq.size(); ++i) {
            std::cout << xq[i] << "\t\t" << yq[i] << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Exception occurred: " << e.what() << std::endl;
    }

    return 0;
}