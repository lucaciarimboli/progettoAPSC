#include <iostream>
#include <vector>
#include "BulkData.hpp"
#include "MeanData.hpp"

int main() {
    // Create some sample data for testing
    std::vector<double> time_data = {0.0, 1.0, 2.0, 3.0, 4.0};
    unsigned int count_sst = 3;
    double N = 100.0;

    // Create some sample MeanData objects
    std::vector<MeanData> mean_data;
    for (int i = 0; i < 5; ++i) {
        MeanData md;
        md.set_position({static_cast<double>(i), static_cast<double>(i + 1), static_cast<double>(i + 2)});
        md.set_sigma({static_cast<double>(i * 0.1), static_cast<double>((i + 1) * 0.1), static_cast<double>((i + 2) * 0.1)});
        mean_data.push_back(md);
    }

    // Create a BulkData object and update it with the sample data
    BulkData bulk_data;
    bulk_data.update_bulk(time_data, count_sst, mean_data, N);

    // Output the results
    std::array<double, 3> w = bulk_data.get_w();
    std::array<double, 3> w_err = bulk_data.get_w_err();
    std::array<double, 3> DN = bulk_data.get_DN();
    std::array<double, 3> DN_err = bulk_data.get_DN_err();

    std::cout << "Drift velocity (w): ";
    for (const auto& val : w) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    std::cout << "Error in drift velocity (w_err): ";
    for (const auto& val : w_err) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    std::cout << "Diffusion constant (DN): ";
    for (const auto& val : DN) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    std::cout << "Error in diffusion constant (DN_err): ";
    for (const auto& val : DN_err) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    return 0;
}