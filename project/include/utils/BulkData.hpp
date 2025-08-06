#ifndef BULKDATA_H
#define BULKDATA_H

#include <array>
#include <vector>
#include <algorithm>

#include "utils/MeanData.hpp"

class BulkData {
    public:

    // Constructor
    BulkData() : w({0.0, 0.0, 0.0}), w_err({0.0, 0.0, 0.0}),
                 DN({0.0, 0.0, 0.0}), DN_err({0.0, 0.0, 0.0}) {};  

    // Public methods:
    void update_bulk(const std::vector<double> & tt, const unsigned int & count_sst, const std::vector<MeanData> & mea, const double & N);
    bool is_empty() const;

    // Getters:
    const std::array<double, 3> & get_w() const { return w; };
    const std::array<double, 3> & get_w_err() const { return w_err; };  
    const std::array<double, 3> & get_DN() const { return DN; };
    const std::array<double, 3> & get_DN_err() const { return DN_err; };

    // Setters:
    void set_mean(const std::vector<MeanData> & y){mean = y;};

    private:

    // Class Members
    std::array<double, 3> w;      // Drift velocity
    std::array<double, 3> w_err;  // Error in drift velocity
    std::array<double, 3> DN;     // Diffusion constant
    std::array<double, 3> DN_err; // Error in diffusion constant
                 
    std::vector<double> t;        // time starting from T_sst
    double t_max;                 // maximum time after steady state
    std::vector<MeanData> mean;   // mean data of steady states

    // Private Methods:
    void normalize(std::vector<double> & y, const double & y_max);
    void update_time_vector(const std::vector<double> & tt, const unsigned int & count_sst);
    void update_mean_data(const unsigned int & count_sst, const std::vector<MeanData> & mea);
    void compute_drift_velocity();
    void compute_diffusion_coeff(const double& N);
    const std::array<double,2> linear_regression(const std::vector<double>& y) const;
};

#endif