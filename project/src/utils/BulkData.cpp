#include "utils/BulkData.hpp"

void BulkData::normalize(std::vector<double> & y, const double & y_max) {
        std::transform(y.begin(), y.end(), y.begin(), [y_max](double yy) { return yy / y_max; });
    };

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
    t_max = t.back();
    normalize(t, t_max);
}

void BulkData::update_mean_data(const unsigned int & count_sst, const std::vector<MeanData> & mm) {
    mean.assign(mm.cend() - count_sst, mm.cend());
    const std::array<double, 3>& pos0 = mean[0].get_position();
    const std::array<double, 3>& var0 = mean[0].get_variance();

    for (auto it = mean.begin() + 1; it != mean.end(); it++) {
        std::array<double, 3> pos = it->get_position();
        std::array<double, 3> var = it->get_variance();

        // Center data around the mean data point at steady state time.
        pos[0] -= pos0[0];
        var[0] -= var0[0];

        pos[1] -= pos0[1];
        var[1] -= var0[1];

        pos[2] -= pos0[2];
        var[2] -= var0[2];

        it->set_position(pos);
        it->set_variance(var);
    }

    mean[0].set_position({0.0, 0.0, 0.0});
    mean[0].set_variance({0.0, 0.0, 0.0});
}

void BulkData::compute_drift(){
    // Each component of y is a vector with the corresponding component of the position for each time after sst.
    std::vector<double> x, y, z;
    x.reserve(mean.size());
    y.reserve(mean.size());
    z.reserve(mean.size());

    for (auto it = mean.cbegin(); it != mean.cend(); it++) {
        const std::array<double, 3>& pos = it->get_position();
        x.push_back(pos[0]);
        y.push_back(pos[1]);
        z.push_back(pos[2]);
    }

    const double x_max = *std::max_element(x.cbegin(), x.cend());
    const double y_max = *std::max_element(y.cbegin(), y.cend());
    const double z_max = *std::max_element(z.cbegin(), z.cend());

    normalize(x, x_max);
    normalize(y, y_max);
    normalize(z, z_max);

    const std::array<double, 2> m0 = linear_regression(x);
    const std::array<double, 2> m1 = linear_regression(y);
    const std::array<double, 2> m2 = linear_regression(z);

    w[0] = m0[0] * x_max / t_max;
    w_err[0] = 0.25 * m0[1] * y_max / t_max;

    w[1] = m1[0] * y_max / t_max;
    w_err[1] = 0.25 * m1[1] * y_max / t_max;

    w[2] = m2[0] * z_max / t_max;
    w_err[2] = 0.25 * m2[1] * z_max / t_max;

}

void BulkData::compute_diffusion(const double N){

    // Each component of y is a vector with the corresponding component of the position for each time after sst.
    std::vector<double> x, y, z;
    x.reserve(mean.size());
    y.reserve(mean.size());
    z.reserve(mean.size());

    for (auto it = mean.cbegin(); it != mean.cend(); it++) {
        const std::array<double, 3>& var = it->get_variance();
        x.push_back(var[0] / 2);
        y.push_back(var[1] / 2);
        z.push_back(var[2] / 2);
    }

    const double x_max = *std::max_element(x.cbegin(), x.cend());
    const double y_max = *std::max_element(y.cbegin(), y.cend());
    const double z_max = *std::max_element(z.cbegin(), z.cend());

    normalize(x, x_max);
    normalize(y, y_max);
    normalize(z, z_max);

    const std::array<double, 2> m0 = linear_regression(x);
    const std::array<double, 2> m1 = linear_regression(y);
    const std::array<double, 2> m2 = linear_regression(z);

    DN[0] = N * m0[0] * x_max / t_max;
    DN_err[0] = N * 0.25 * m0[1] * y_max / t_max;

    DN[1] = N * m1[0] * y_max / t_max;
    DN_err[1] = N * 0.25 * m1[1] * y_max / t_max;

    DN[2] = N * m2[0] * z_max / t_max;
    DN_err[2] = N * 0.25 * m2[1] * z_max / t_max;
}

const std::array<double,2> BulkData::linear_regression(const std::vector<double>& y) const{
        
    size_t n = y.size();

    // Compute mean values:
    const double x_mean = std::accumulate(t.cbegin(), t.cend(), 0.0) / n;
    const double y_mean = std::accumulate(y.cbegin(), y.cend(), 0.0) / n;

    // Compute covariance and variance:
    std::vector<double> xx(n);
    std::vector<double> yy(n);
    std::transform(t.cbegin(), t.cend(), xx.begin(), [x_mean](double x) { return x - x_mean; });
    std::transform(y.cbegin(), y.cend(), yy.begin(), [y_mean](double y) { return y - y_mean; });
    const double cov = std::inner_product(xx.cbegin(), xx.cend(), yy.cbegin(), 0.0);
    const double var = std::inner_product(xx.cbegin(), xx.cend(), xx.cbegin(), 0.0);

    // Compute slope and intercept:
    const double m = cov / var;
    const double q = y_mean - m * x_mean;

    // Compute residuals:
    std::vector<double> res(n);
    std::transform(t.cbegin(), t.cend(), y.cbegin(), res.begin(),
        [m, q](double x, double y) { return y - (q + m * x); }
    );

    // Compute standard error for slope:
    const double standard_err= std::sqrt(std::inner_product(res.cbegin(), res.cend(), res.cbegin(), 0.0) / ((n - 2) * var));

    // Confidence interval of 95% for slope:

    // The smallest # of dof is n = 8. t_value = 1.860 for 95% confidence interval with 8 dof.
    // A t-student quantile table would be required, to avoid importing external libraries just for this task,
    // an estimate using a linear behavior of the quantile from 8 dof (1.860) to infinite dof (1.645) is applied.

    // Consider that this is a rough estimate but it is accetable since the convergence is expected to
    // verify roughly at n of the order of 10^5, so the actual "t_value" becomes after few iterations
    // very close to 1.645.

    // This estimate is safe as the actual t-student values decrease faster than linearly
    // (e.g. for 100 d.o.f: actual t_value = 1.660, computed value: 1.662)

    // Confidence interval of 95% for slope with 8 dof (minimum value of n-2):
    //double t_value = 1.860; 
    double t_value = 1.645 + 1.72 / (n - 2); // linear interpolation for t-student quantile
    return { m, 2 * t_value * standard_err };
};
