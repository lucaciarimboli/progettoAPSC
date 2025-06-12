#include "utils/FluxData.hpp"

FluxData::FluxData():
    v_int_sum({0.0, 0.0, 0.0}),
    w({0.0, 0.0, 0.0}),
    D_sum({0.0, 0.0, 0.0}),
    DN({0.0, 0.0, 0.0})
{};

void FluxData::compute_drift_velocity(const mc::MATRIX & v_int, const double t_total){
    for(auto it = v_int.cbegin(); it != v_int.cend(); it++){
        v_int_sum[0] += (*it)[0];
        v_int_sum[1] += (*it)[1];
        v_int_sum[2] += (*it)[2];   
    }

    w[0] = v_int_sum[0] / t_total;
    w[1] = v_int_sum[1] / t_total;
    w[2] = v_int_sum[2] / t_total;
};

void FluxData::compute_diffusion_const(const mc::MATRIX & r, const mc::MATRIX & v,
    const double den, const double count_sst){

    // Compute mean (in time) position, velocity and diffusion flux:
    const std::array<double,3> mean_r = compute_mean(r);
    const std::array<double,3> mean_v = compute_mean(v);
    const std::array<double,3> D_flux = compute_mean(elementwise_product(r,v));

    // Sum diffusion constants over all iterations:
    D_sum[0] += D_flux[0] - mean_r[0] * mean_v[0];
    D_sum[1] += D_flux[1] - mean_r[1] * mean_v[1];
    D_sum[2] += D_flux[2] - mean_r[2] * mean_v[2];

    // Update diffusion constant:
    const size_t N = count_sst - 10;    // Number of iterations since flux data are updated
    DN[0] = den * D_sum[0] / N;
    DN[1] = den * D_sum[1] / N;
    DN[2] = den * D_sum[2] / N;
};

const std::array<double,3> FluxData::compute_mean(const mc::MATRIX & M) const{
    // Compute component-by-component mean
    std::array<double,3> mean({0.0, 0.0, 0.0});
    const size_t n = M.size();

    for(auto it = M.cbegin(); it != M.cend(); it++){
        mean[0] += (*it)[0];
        mean[1] += (*it)[1];
        mean[2] += (*it)[2];   
    }

    mean[0] /= n;
    mean[1] /= n;
    mean[2] /= n;

    return mean;
};

const mc::MATRIX FluxData::elementwise_product(const mc::MATRIX & A, const mc::MATRIX & B) const {

    // Equivalent of Matlab's .* product between matrices
    // Remark that MATRIX = std::vector<std::array<double,3>>
    const size_t n = A.size();

    // Size comparison for safety reasons
    //if( n != B.size()){
    //    throw std::invalid_argument("Matrices must have the same size");
    //}

    // Compute product element-by-element:
    mc::MATRIX A_dot_B;
    A_dot_B.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        A_dot_B.emplace_back(std::array<double, 3>{
            A[i][0] * B[i][0],
            A[i][1] * B[i][1],
            A[i][2] * B[i][2]
        });
    }

    return A_dot_B;
};