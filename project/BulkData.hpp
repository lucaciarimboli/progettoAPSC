#ifndef BULKDATA_H
#define BULKDATA_H

#include <array>
#include <vector>
#include <algorithm>
#include "MeanData.hpp"

class MeanData;

class BulkData {
    public:
    // Constructor
    BulkData() : w({0.0, 0.0, 0.0}), w_err({0.0, 0.0, 0.0}),
                 DN({0.0, 0.0, 0.0}), DN_err({0.0, 0.0, 0.0}) {};
    
    // PER MIGLIORARE update_bulk() CONSIDERA DI AGGIORNARE SOLO L'ELEMENTO IN CODA DI t E mean --> DA RIVEDERE ALLA FINE
    void update_bulk(const std::vector<double> & tt, const unsigned int & count_sst, const std::vector<MeanData> & mea, const double & N){
        update_time_vector(tt, count_sst);
        update_mean_data(count_sst, mea);
        compute_drift();
        compute_diffusion(N);
    };

    // Getters:
    const std::array<double, 3> & get_w() const { return w; };
    const std::array<double, 3> & get_w_err() const { return w_err; };  
    const std::array<double, 3> & get_DN() const { return DN; };
    const std::array<double, 3> & get_DN_err() const { return DN_err; };
    // Setters:
    void set_t(const std::vector<double> & tt){
         t = tt;
         t_max = *std::max_element(t.cbegin(), t.cend());
         normalize(t,t_max);
    };
    void set_mean(const std::vector<MeanData> & y){mean = y;};

    private:
    std::array<double, 3> w;      // Drift velocity
    std::array<double, 3> w_err;  // Error in drift velocity
    std::array<double, 3> DN;     // Diffusion constant
    std::array<double, 3> DN_err; // Error in diffusion constant

    double t_max;                 // maximum of t --> SEMPLICEMENTE t.end() ??
    std::vector<double> t;        // time starting from T_sst
    std::vector<MeanData> mean;   // mean data of steady states

    void normalize(std::vector<double> & y, const double & y_max) {
        std::transform(y.begin(), y.end(), y.begin(), [y_max](double yy) { return yy / y_max; });
    };

    // Multiple variables linear regression, returns the coefficients of the linear regression and uncertainty.
    // Returned array is B = [m,u_m] where:
    // y[i] = q + m*t[i] + res, m € (m-um/2,m+um/2) with >95% confidence.
    std::array<double,2> linear_regression(const std::vector<double>& y){
        
        size_t n = y.size();

        // Compute mean values:
        double x_mean = std::accumulate(t.cbegin(), t.cend(), 0.0) / n;
        double y_mean = std::accumulate(y.cbegin(), y.cend(), 0.0) / n;

        // Compute covariance and variance:
        std::vector<double> xx(n);
        std::vector<double> yy(n);
        std::transform(t.cbegin(), t.cend(), xx.begin(), [x_mean](double x) { return x - x_mean; });
        std::transform(y.cbegin(), y.cend(), yy.begin(), [y_mean](double y) { return y - y_mean; });
        double cov = std::inner_product(xx.cbegin(), xx.cend(), yy.cbegin(), 0.0);
        double var = std::inner_product(xx.cbegin(), xx.cend(), xx.cbegin(), 0.0);

        // Compute slope and intercept:
        double m = cov / var;
        double q = y_mean - m * x_mean;

        // Compute residuals:
        std::vector<double> res(n);
        std::transform(t.cbegin(), t.cend(), y.cbegin(), res.begin(), [m, q](double x, double y) { return y - (q + m * x); });

        // Compute standard error for slope:
        double standard_err= std::sqrt(std::inner_product(res.cbegin(), res.cend(), res.cbegin(), 0.0) / ((n - 2) * var));

        // Confidence interval of 95% for slope:

        // The smallest # of dof is n = 8. t_value = 1.860 for 95% confidence interval with 8 dof.
        // A t-student quantile table would be required, to avoid importing external libraries just for this task,
        // an overestimate of the confidence interval is acceptable.

        // Confidence interval of 95% for slope with 8 dof (minimum value of n-2):
        double t_value = 1.860; 
        return { m, 2 * t_value * standard_err };
    };

    void update_time_vector(const std::vector<double> & tt, const unsigned int & count_sst) {
        t.assign(tt.cend() - count_sst, tt.cend());
        std::transform(t.begin() + 1, t.end(), t.begin() + 1, [this](double val) { return val - t[0]; });
        t[0] = 0.0;
        t_max = *std::max_element(t.cbegin(), t.cend());
        normalize(t, t_max);
    }

    void update_mean_data(const unsigned int & count_sst, const std::vector<MeanData> & mea) {
        mean.assign(mea.cend() - count_sst, mea.cend());
        const std::array<double, 3> & pos0 = mean[0].get_position();
        const std::array<double, 3> & sig0 = mean[0].get_sigma();

        for (auto it = mean.begin() + 1; it != mean.end(); it++) {
            std::array<double, 3> pos = it->get_position();
            std::array<double, 3> sig = it->get_sigma();

            for (size_t j = 0; j < 3; j++) {
                pos[j] -= pos0[j];
                sig[j] = 0.5 * std::pow(sig[j], 2) - 0.5 * std::pow(sig0[j], 2);
            }
        }

        mean[0].set_position({0.0, 0.0, 0.0});
        mean[0].set_sigma({0.0, 0.0, 0.0});
    }

    void compute_drift(){
        // Each component of y is a vector with the corresponding component of the position for each time after sst.
        std::array<std::vector<double>, 3> y;
        for (const MeanData& m : mean) {
            std::array<double, 3> pos = m.get_position();
            for (size_t i = 0; i < 3; i++) {
                y[i].emplace_back(pos[i]);
            }
        }

        for (size_t i = 0; i < 3; i++) {
            double y_max = *std::max_element(y[i].cbegin(), y[i].cend());
            normalize(y[i], y_max);
            std::array<double, 2> m = linear_regression(y[i]);
            w[i] = m[0] * y_max / t_max;
            w_err[i] = 0.25 * m[1] * y_max / t_max;
        }
    }

    void compute_diffusion(const double & N){
        // Each component of y is a vector with the corresponding component of the position for each time after sst.
        std::array<std::vector<double>, 3> y;
        for (const MeanData& m : mean) {
            std::array<double, 3> sigma = m.get_sigma();
            for (size_t i = 0; i < 3; i++) {
                y[i].emplace_back(sigma[i]);
            }
        }

        for (size_t i = 0; i < 3; i++) {
            double y_max = *std::max_element(y[i].cbegin(), y[i].cend());
            normalize(y[i], y_max);
            std::array<double, 2> m = linear_regression(y[i]);
            DN[i] = N * m[0] * y_max / t_max;
            DN_err[i] = N * 0.25 * m[1] * y_max / t_max;
        }
    }
};

#endif