#ifndef BULKDATA_H
#define BULKDATA_H

#include <array>
#include <vector>
#include <algorithm>

#include "utils/MeanData.hpp"

/**
 * @class BulkData
 * @brief Computes and stores bulk transport coefficients of the electron swarm.
 *
 * This class manages the computation of the bulk drift velocity and diffusion coefficient
 * from stady-state mean data obtained in a Monte Carlo simulation.
 * The computation is performed by linear regression.
 * The errors on bulk velocity and diffusion coefficent are computed through the 
 * 95% confidence interval on the slope coefficient from the linear regression.
 */
class BulkData {
    public:

    /**
     * @brief Constructor.
     */
    BulkData() : w({0.0, 0.0, 0.0}), w_err({0.0, 0.0, 0.0}),
                 DN({0.0, 0.0, 0.0}), DN_err({0.0, 0.0, 0.0}) {};  

    // Public methods:

    /**
     * @brief Updates the bulk transport coefficients.
     * 
     * Uses the time vector and averaged data after the steady state to compute
     * drift velocity, diffusion coefficient and corresponding errors.
     *
     * @param tt Vector of time steps.
     * @param count_sst Number of steady state time steps.
     * @param mea Mean data over time.
     * @param N Gas number density [m^-3].
     */
    void update_bulk(const std::vector<double> & tt, const unsigned int & count_sst, const std::vector<MeanData> & mea, const double & N);

    /**
     * @brief Checks if no bulk data has been computed yet.
     * 
     * @return true if all arrays are zeros, false otherwise.
     */
    bool is_empty() const;

    // Getters:

    /**
     * @brief Returns the bulk drift velocity vector.
     * @return Array with drift velocity [m/s].
     */
    const std::array<double, 3> & get_w() const { return w; };

    /**
     * @brief Returns the drift velocity error.
     * @return Array with drift velocity uncertainty [m/s].
     */
    const std::array<double, 3> & get_w_err() const { return w_err; };

    /**
     * @brief Returns the diffusion coefficient.
     * @return Array with diffusion coefficient [m²/s].
     */
    const std::array<double, 3> & get_DN() const { return DN; };

    /**
     * @brief Returns the error on the diffusion constant.
     * @return Array with diffusion coefficient uncertainty [m²/s].
     */
    const std::array<double, 3> & get_DN_err() const { return DN_err; };

    // Setters:

    /**
     * @brief Stores the time-averaged MeanData after equilibrium.
     *
     * @param y Vector of MeanData values to store.
     */
    void set_mean(const std::vector<MeanData> & y){mean = y;};

    private:

    // Class Members
    std::array<double, 3> w;      ///< Drift velocity
    std::array<double, 3> w_err;  ///< Error in drift velocity
    std::array<double, 3> DN;     ///< Diffusion constant
    std::array<double, 3> DN_err; ///< Error in diffusion constant
                 
    std::vector<double> t;        ///< time starting from T_sst
    double t_max;                 ///< maximum time after steady state
    std::vector<MeanData> mean;   ///< mean data of steady states

    // Private Methods:

    /**
     * @brief Normalizes the input vector between 0 and 1 by y_max.
     *
     * @param y Vector to normalize.
     * @param y_max Maximum value for normalization.
     */
    void normalize(std::vector<double> & y, const double & y_max);

    /**
     * @brief Extracts times after steady state.
     *
     * @param tt Sull simulation time vector.
     * @param count_sst Number of steady state time steps.
     */
    void update_time_vector(const std::vector<double> & tt, const unsigned int & count_sst);

    /**
     * @brief Extracts mean data of steady state time steps.
     *
     * @param count_sst Number of time steps at steady state.
     * @param mea Vector of MeanData from the simulation.
     */
    void update_mean_data(const unsigned int & count_sst, const std::vector<MeanData> & mea);

    /**
     * @brief Computes the drift velocity using regression on electrons center of mass.
     * 
     * The center of mass derivative is approximated using the slope of the linear regression
     * of the center of mass against the time.
     * The linear regression considers only steady state data, both the center of mass and
     * the time vector are normalized to 1.
     */
    void compute_drift_velocity();

    /**
     * @brief Computes the drift velocity using regression on electrons center of mass.
     * @param N Gas number density [m^-3].
     * 
     * The variance derivative is approximated using the slope of the linear regression
     * of the  variance in position against the time.
     * The linear regression considers only steady state data, both the variance in position and
     * the time vector are normalized to 1.
     */
    void compute_diffusion_coeff(const double& N);

    /**
     * @brief Performs linear regression against a normalized time vector.
     *
     * @param y Vector of data values.
     * @return Array with slope and 95% confidence interval width on slope.
     */
    const std::array<double,2> linear_regression(const std::vector<double>& y) const;
};

#endif