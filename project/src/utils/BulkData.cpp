#include "utils/BulkData.hpp"

void BulkData::update_bulk(const std::vector<double> & tt, const unsigned int & count_sst, const std::vector<MeanData> & mea, const double & N){
    update_time_vector(tt, count_sst);
    update_mean_data(count_sst, mea);
    compute_drift();
    compute_diffusion(N);
};

void BulkData::update_time_vector(const std::vector<double> & tt, const unsigned int & count_sst) {
    t.assign(tt.cend() - count_sst, tt.cend());
    std::transform(t.begin() + 1, t.end(), t.begin() + 1, [this](double val) { return val - t[0]; });
    t[0] = 0.0;
    normalize(t, t.back());
}

void BulkData::update_mean_data(const unsigned int & count_sst, const std::vector<MeanData> & mea) {
    mean.assign(mea.cend() - count_sst, mea.cend());
    const std::array<double, 3>& pos0 = mean[0].get_position();
    const std::array<double, 3>& sig0 = mean[0].get_sigma();

    for (auto it = mean.begin() + 1; it != mean.end(); it++) {
        std::array<double, 3> pos = it->get_position();
        std::array<double, 3> sig = it->get_sigma();

        for (size_t j = 0; j < 3; j++) {
            pos[j] -= pos0[j];
            sig[j] = 0.5 * (sig[j] * sig[j] - sig0[j] * sig0[j]);;
        }

        it->set_position(pos);
        it->set_sigma(sig);
    }

    mean[0].set_position({0.0, 0.0, 0.0});
    mean[0].set_sigma({0.0, 0.0, 0.0});
}

void BulkData::compute_drift(){
    // Each component of y is a vector with the corresponding component of the position for each time after sst.
    std::array<std::vector<double>, 3> y;
    for (size_t i = 0; i < 3; i++) {
        y[i].reserve(mean.size());
    }
    
    for (const MeanData& m : mean) {
        const std::array<double, 3>& pos = m.get_position();
        for (size_t i = 0; i < 3; i++) {
            y[i].push_back(pos[i]);
        }
    }

    const double t_max = t.back();
    for (size_t i = 0; i < 3; i++) {
        const double y_max = *std::max_element(y[i].cbegin(), y[i].cend());
        normalize(y[i], y_max);
        const std::array<double, 2> m = linear_regression(y[i]);
        w[i] = m[0] * y_max / t_max;
        w_err[i] = 0.25 * m[1] * y_max / t_max;
    }
}

void BulkData::compute_diffusion(const double N){

    // Each component of y is a vector with the corresponding component of the position for each time after sst.
    std::array<std::vector<double>, 3> y;
    for (const MeanData& m : mean) {
        std::array<double, 3> sigma = m.get_sigma();
        for (size_t i = 0; i < 3; i++) {
            y[i].emplace_back(sigma[i]);
        }
    }

    const double t_max = t.back();
    for (size_t i = 0; i < 3; i++) {
        const double y_max = *std::max_element(y[i].cbegin(), y[i].cend());
        normalize(y[i], y_max);
        const std::array<double, 2> m = linear_regression(y[i]);
        DN[i] = N * m[0] * y_max / t_max;
        DN_err[i] = N * 0.25 * m[1] * y_max / t_max;
    }
}

std::array<double,2> BulkData::linear_regression(const std::vector<double>& y){
        
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