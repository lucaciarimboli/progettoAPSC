#ifndef FLUXDATA_H
#define FLUXDATA_H

#include <vector>
#include <array>
#include <cmath>
#include <numeric>

class FluxData {

    public:
    // Constructor:
    FluxData()
        : v_int_sum({0.0, 0.0, 0.0}),
          w({0.0, 0.0, 0.0}),
          D_sum({0.0, 0.0, 0.0}),
          DN({0.0, 0.0, 0.0}),
          N(0) {};

    void compute_drift_velocity(const MATRIX & v_int, const double & t_total){
        for(const std::array<double,3> &v : v_int){
            for( size_t i = 0; i < 3; i++){
                v_int_sum[i] += v[i];
            }
        }

        for( size_t i = 0; i < 3; i++){
            w[i] = v_int_sum[i] / t_total;
        } 
    };

    void compute_diffusion_const(const MATRIX & r, const MATRIX & v, const double & den){
        N++;
        std::array<double,3> mean_r = compute_mean(r);
        std::array<double,3> mean_v = compute_mean(v);
        std::array<double,3> D_flux = compute_mean(dotstar(r,v));

        for(size_t i = 0; i < 3; i++){
            double D_flux_i = D_flux[i] - mean_r[i] * mean_v[i];
            D_sum[i] += D_flux_i;
            DN[i] = den * D_sum[i] / N;     // den is gas density [m^-3]
        }
    }

    private:
    std::array<double, 3> v_int_sum;   // Integrated velocity sum
    std::array<double, 3> w;          // Drift velocity
    std::array<double, 3> D_sum;      // Sum of diffusion constants
    std::array<double, 3> DN;         // Diffusion constant
    size_t N;                         // Number of time steps

    // Compute component-by-component mean
    std::array<double,3> compute_mean(const MATRIX & M){

        // Matrix cannot be empty
        if (M.empty()) {
           throw std::invalid_argument("Matrices cannot be empty.");
        }

        std::array<double,3> mean {};

        for(const std::array<double,3> &v : M){
            for( size_t i = 0; i < 3; i++){
                mean[i] += v[i];
            }
        }

        for( size_t i = 0; i < 3; i++){
            mean[i] = mean[i] / M.size();
        }

        return mean;
    };

    // Equivalent of Matlab's .* product between matrices
    MATRIX dotstar(const MATRIX & A, const MATRIX & B){

        // Matrices cannot be empty
        if (A.empty() || B.empty()) {
           throw std::invalid_argument("Matrices cannot be empty.");
        }

        // Remark that MATRIX = std::<vector<std::array<double,3>>
        size_t n = A.size();

        // Size comparison for safety reasons
        if( n != B.size()){
            throw std::invalid_argument("Matrices must have the same size");
        }

        // Compute product element-by-element:
        MATRIX A_dot_B;
        A_dot_B.reserve(n);
        for( size_t i = 0; i < n; i++){
            std::array<double,3> column;
            for( size_t j = 0; j < 3; j++){
                column[j] = A[i][j] * B[i][j];
            }
            A_dot_B.emplace_back(column);
        }

        return A_dot_B;
    };

};

#endif