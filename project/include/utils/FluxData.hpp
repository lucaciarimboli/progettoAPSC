#ifndef FLUXDATA_H
#define FLUXDATA_H

#include <vector>
#include <array>
#include <cmath>
#include <numeric>

#include "Common.hpp"

class FluxData {

    public:
    // Constructor:
    FluxData();

    // Public methods:
    void compute_drift_velocity(const mc::MATRIX& v_int, const double& t_total);
    void compute_diffusion_const(const mc::MATRIX& r, const mc::MATRIX& v, const double& den, const double& count_sst);

    // Setters:
    void set_drift_vel(const std::array<double,3>& ww){ w = ww; };
    void set_D(const std::array<double,3>& D){ DN = D; };

    // Getters:
    const std::array<double,3> & get_w() const { return w; };
    const std::array<double,3> & get_DN() const { return DN; };

    private:

    // Class members:
    std::array<double, 3> v_int_sum;  // Integrated velocity sum
    std::array<double, 3> w;          // Drift velocity (z-direction)

    std::array<double, 3> D_sum;      // Sum of diffusion constants
    std::array<double, 3> DN;         // Diffusion constant

    // Private methods:
    const std::array<double,3> compute_mean(const mc::MATRIX& M) const;
    const mc::MATRIX elementwise_product(const mc::MATRIX& A, const mc::MATRIX& B) const;
};

#endif