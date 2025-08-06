#ifndef FLUXDATA_H
#define FLUXDATA_H

#include <vector>
#include <array>
#include <cmath>
#include <numeric>

#include "Common.hpp"

/**
 * @class FluxData
 * @brief Computes and stores flux transport coefficients of the electron swarm.
 *
 * This class manages the calculation of macroscopic electron transport parameters such as 
 * drift velocity and diffusion coefficient using integrated microscopic data over time.
 */
class FluxData {

    public:
    
    /**
     * @brief Constructor.
     */
    FluxData();

    // Public methods:

    /**
     * @brief Computes the drift velocity vector.
     * 
     * @param v_int Integrated velocity vector over time for each electron.
     * @param t_total Total integration time.
     */
    void compute_drift_velocity(const mc::MATRIX& v_int, const double& t_total);

    /**
     * @brief Computes the diffusion constant in each direction.
     * 
     * @param r Electron positions over time.
     * @param v Electron velocities over time.
     * @param den Background gas number density.
     * @param count_sst Number of steady state time steps.
     */
    void compute_diffusion_const(const mc::MATRIX& r, const mc::MATRIX& v, const double& den, const double& count_sst);

    // Getters:

    /**
     * @brief Returns the drift velocity vector.
     * 
     * @return Drift velocity vector.
     */
    const std::array<double,3> & get_w() const { return w; };
    /**
     * @brief Returns the diffusion coefficients vector.
     * 
     * @return Diffusion coefficients vector.
     */
    const std::array<double,3> & get_DN() const { return DN; };

    private:

    // Class members:
    std::array<double, 3> v_int_sum;  ///< Integrated velocity sum
    std::array<double, 3> w;          ///< Drift velocity (z-direction)

    std::array<double, 3> D_sum;      ///< Sum of diffusion constants
    std::array<double, 3> DN;         ///< Diffusion constant

    // Private methods:
    
    /**
     * @brief Computes the mean of a given matrix along its first index.
     * 
     * @param M The matrix whose column-wise mean is to be computed.
     * @return The resulting 3D mean vector.
     */
    const std::array<double,3> compute_mean(const mc::MATRIX& M) const;

    /**
     * @brief Performs element-wise product between two matrices of the same size.
     * 
     * @param A First matrix.
     * @param B Second matrix.
     * @return The resulting matrix A .* B .
     */
    const mc::MATRIX elementwise_product(const mc::MATRIX& A, const mc::MATRIX& B) const;
};

#endif