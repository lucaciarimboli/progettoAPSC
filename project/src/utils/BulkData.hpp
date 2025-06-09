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
    
    // PER MIGLIORARE update_bulk() CONSIDERA DI AGGIORNARE SOLO L'ELEMENTO IN CODA DI t E mean --> DA RIVEDERE ALLA FINE
    void update_bulk(const std::vector<double> & tt, const unsigned int & count_sst, const std::vector<MeanData> & mea, const double & N);

    // Getters:
    const std::array<double, 3> & get_w() const { return w; };
    const std::array<double, 3> & get_w_err() const { return w_err; };  
    const std::array<double, 3> & get_DN() const { return DN; };
    const std::array<double, 3> & get_DN_err() const { return DN_err; };
    // Setters:
    /*void set_t(const std::vector<double> & tt){
         t = tt;
         t_max = t.back();
         normalize(t,t_max);
    };*/
    void set_mean(const std::vector<MeanData> & y){mean = y;};

    // Check if bulk has meaningful data or not yet (steady state not reached or not enough data):
    bool is_empty() const {
        return (w[0] == 0.0 && w[1] == 0.0 && w[2] == 0.0 &&
                DN[0] == 0.0 && DN[1] == 0.0 && DN[2] == 0.0 &&
                w_err[0] == 0.0 && w_err[1] == 0.0 && w_err[2] == 0.0 &&
                DN_err[0] == 0.0 && DN_err[1] == 0.0 && DN_err[2] == 0.0);
    };

    private:
    std::array<double, 3> w;      // Drift velocity
    std::array<double, 3> w_err;  // Error in drift velocity
    std::array<double, 3> DN;     // Diffusion constant
    std::array<double, 3> DN_err; // Error in diffusion constant
                 
    std::vector<double> t;        // time starting from T_sst
    double t_max;                 // maximum time after steady state
    std::vector<MeanData> mean;   // mean data of steady states

    void normalize(std::vector<double> & y, const double & y_max) {
        std::transform(y.begin(), y.end(), y.begin(), [y_max](double yy) { return yy / y_max; });
    };

    void update_time_vector(const std::vector<double> & tt, const unsigned int & count_sst);
    void update_mean_data(const unsigned int & count_sst, const std::vector<MeanData> & mea);
    void compute_drift();
    void compute_diffusion(const double N);

    // Multiple variables linear regression, returns the coefficients of the linear regression and uncertainty.
    // Returned array is B = [m,u_m] where:
    // y[i] = q + m*t[i] + res, m € (m-um/2,m+um/2) with >95% confidence.
    std::array<double,2> linear_regression(const std::vector<double>& y);
};

#endif