#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <string>
#include <cmath>    // For std::lerp --> IT REQUIRES C++20 !

// This class can be used to perform interpolation on a set of data points. 
// This class has been developed when needing a 1D linear interpolation, which is so far the only kind of interpolation supported.
// However, the class has been designed to be extended to support other kinds of interpolation in the future.

class Interpolation {
public:

    // Constructor
    Interpolation() = default;
    // Constructor with kind
    Interpolation(const std::string& kind) {
        setInterpolationKind(kind);
    }

    // Interpolation function
    std::vector<double> interpolation(const std::vector<double>& x,
                                 auto& y,
                                 const std::vector<double>& xq) {
        if (interpolation_kind == "linear") {   // Add other kinds as needed
            return linear_interpolation(x, y, xq);
        } else {
            throw std::invalid_argument("Unsupported interpolation kind");
        }
    }

    // Set the interpolation kind
    void setInterpolationKind(const std::string& kind) {
        if (kind == "linear") {     // Add other kinds as needed
            interpolation_kind = kind;
        } else {
            throw std::invalid_argument("Unsupported interpolation kind");
        }
    }
    // Get the current interpolation kind
    std::string getInterpolationKind() const {
        return interpolation_kind;
    }


private:
    std::string interpolation_kind = "linear"; // Default interpolation kind

    // Such class might be expanded to support other interpolation methods
    // like cubic, spline, etc. in the future.
    // Multi-column interpolation: each row of v corresponds to x, each column a dataset
    static std::vector<double> linear_interpolation(const std::vector<double>& x,
                                       const std::vector<std::vector<double>>& y,
                                       const std::vector<double>& xq) {
        if (x.size() != y.size()) {
            throw std::invalid_argument("x and v must have the same number of rows");
        }

        size_t num_points = xq.size();
        size_t num_cols = y[0].size();
        std::vector<double> result(num_points * num_cols);

        for (size_t n = 0; n < num_points; n++) {

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

            for (size_t col = 0; col < num_cols; ++col) {
                double yi0 = y[i][col];
                double yi1 = y[i + 1][col];
                result[n * num_cols + col] = std::lerp(yi0, yi1, t);
            }
        }

        return result;
    }

    // Single column interpolation
    static std::vector<double> linear_interpolation(const std::vector<double>& x,
                                       const std::vector<double>& y,
                                       const std::vector<double>& xq) {
        std::vector<std::vector<double>> matrix_y(y.size(), std::vector<double>(1));
        for (size_t i = 0; i < y.size(); i++)
            matrix_y[i][0] = y[i];

        std::vector<double> flat_result = linear_interpolation(x, matrix_y, xq);

        std::vector<double> result(xq.size());
        for (size_t i = 0; i < xq.size(); i++)
            result[i] = flat_result[i];

        return result;
    }
};

#endif // INTERPOLATION_HPP