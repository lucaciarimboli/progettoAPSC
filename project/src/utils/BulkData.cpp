#include "utils/BulkData.hpp"

void BulkData::normalize(std::vector<double> & y, const double & y_max) {
    std::transform(y.begin(), y.end(), y.begin(), [y_max](const double& yy) { return yy / y_max; });
};

void BulkData::update_bulk(const std::vector<double> & tt, const unsigned int & count_sst, const std::vector<MeanData> & mea, const double & N){
    update_time_vector(tt, count_sst);
    update_mean_data(count_sst, mea);
    compute_drift_velocity();
    compute_diffusion_coeff(N);
};

bool BulkData::is_empty() const {
    // Check if bulk has meaningful data or not yet (steady state not reached or not enough data):
    return (w[0] == 0.0 && w[1] == 0.0 && w[2] == 0.0 &&
            DN[0] == 0.0 && DN[1] == 0.0 && DN[2] == 0.0 &&
            w_err[0] == 0.0 && w_err[1] == 0.0 && w_err[2] == 0.0 &&
            DN_err[0] == 0.0 && DN_err[1] == 0.0 && DN_err[2] == 0.0);
};

void BulkData::update_time_vector(const std::vector<double> & tt, const unsigned int & count_sst) {
    t.assign(tt.cend() - count_sst, tt.cend());
    const double t0 = t[0];
    std::transform(t.begin(), t.end(), t.begin(), [t0](const double& val) { return val - t0; });
    t_max = t.back();
    normalize(t, t_max);
}

void BulkData::update_mean_data(const unsigned int & count_sst, const std::vector<MeanData> & mm) {
    mean.assign(mm.cend() - count_sst, mm.cend());
    const std::array<double, 3> pos0 = mean[0].get_position();
    const std::array<double, 3> var0 = mean[0].get_variance();

    for (auto it = mean.begin(); it != mean.end(); it++) {
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
}

void BulkData::compute_drift_velocity(){
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

    const double x_max = *std::max_element(x.cbegin(), x.cend(), 
        [](const double& a, const double& b){return std::less{}(std::abs(a),std::abs(b));}
    );
    const double y_max = *std::max_element(y.cbegin(), y.cend(),
        [](const double& a, const double& b){return std::less{}(std::abs(a),std::abs(b));}
    );
    const double z_max = *std::max_element(z.cbegin(), z.cend(),
        [](const double& a, const double& b){return std::less{}(std::abs(a),std::abs(b));}
    );

    normalize(x, x_max);
    normalize(y, y_max);
    normalize(z, z_max);

    const std::array<double, 2> m0 = linear_regression(x);
    const std::array<double, 2> m1 = linear_regression(y);
    const std::array<double, 2> m2 = linear_regression(z);

    w[0] = m0[0] * x_max / t_max;
    w_err[0] = 0.25 * m0[1] * x_max / t_max;

    w[1] = m1[0] * y_max / t_max;
    w_err[1] = 0.25 * m1[1] * y_max / t_max;

    w[2] = m2[0] * z_max / t_max;
    w_err[2] = 0.25 * m2[1] * z_max / t_max;
}

void BulkData::compute_diffusion_coeff(const double& N){

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

    const double x_max = *std::max_element(x.cbegin(), x.cend(), 
        [](const double& a, const double& b){return std::less{}(std::abs(a),std::abs(b));}
    );
    const double y_max = *std::max_element(y.cbegin(), y.cend(),
        [](const double& a, const double& b){return std::less{}(std::abs(a),std::abs(b));}
    );
    const double z_max = *std::max_element(z.cbegin(), z.cend(),
        [](const double& a, const double& b){return std::less{}(std::abs(a),std::abs(b));}
    );

    normalize(x, x_max);
    normalize(y, y_max);
    normalize(z, z_max);

    const std::array<double, 2> m0 = linear_regression(x);
    const std::array<double, 2> m1 = linear_regression(y);
    const std::array<double, 2> m2 = linear_regression(z);

    DN[0] = N * m0[0] * x_max / t_max;
    DN_err[0] = N * 0.25 * m0[1] * x_max / t_max;

    DN[1] = N * m1[0] * y_max / t_max;
    DN_err[1] = N * 0.25 * m1[1] * y_max / t_max;

    DN[2] = N * m2[0] * z_max / t_max;
    DN_err[2] = N * 0.25 * m2[1] * z_max / t_max;
}

const std::array<double,2> BulkData::linear_regression(const std::vector<double>& y) const{
    // Multiple variables linear regression, returns the coefficients of the linear regression and uncertainty.
    // Returned array is B = [m,IC] where:
    // y[i] = q + m*t[i] + res, m € (m-IC/2,m+IC/2) with >95% confidence.
        
    const size_t n = y.size();
    
    // Compute mean values:
    const double x_mean = std::accumulate(t.cbegin(), t.cend(), 0.0) / n;
    const double y_mean = std::accumulate(y.cbegin(), y.cend(), 0.0) / n;

    // Compute covariance and variance:
    std::vector<double> xx(n);
    std::vector<double> yy(n);
    std::transform(t.cbegin(), t.cend(), xx.begin(), [x_mean](const double& x) { return x - x_mean; });
    std::transform(y.cbegin(), y.cend(), yy.begin(), [y_mean](const double& y) { return y - y_mean; });
    const double cov = std::inner_product(xx.cbegin(), xx.cend(), yy.cbegin(), 0.0);
    const double var = std::inner_product(xx.cbegin(), xx.cend(), xx.cbegin(), 0.0);

    // Compute slope and intercept:
    const double m = cov / var;
    const double q = y_mean - m * x_mean;

    // Compute residuals:
    std::vector<double> res(n);
    std::transform(t.cbegin(), t.cend(), y.cbegin(), res.begin(),
        [m, q](const double& x, const double& y) { return y - (q + m * x); }
    );

    // Compute standard error for slope:
    const double standard_err= std::sqrt(std::inner_product(res.cbegin(), res.cend(), res.cbegin(), 0.0) / ((n - 2) * var));

    // Confidence interval of 95% for slope:

    // The smallest # of dof is n = 9. t_value = 2.2621 for 97.5% confidence interval with 9 dof.
    // A t-student quantile table would be required, to avoid importing external libraries just for this task,
    // an over-estimate using the quantile with 9 dof (2.262) and with infinite dof (1.960) is applied.

    // This estimate is safe as the actual t-student values decrease faster, so this approximation simply provides
    // a slightly more conservative error estimate.
    // (e.g. for 100 d.o.f: actual t_value = 1.984, computed value: 1.9877)

    // Length of confidence interval of 95% for slope with 9 dof (minimum value of n-2):
    const double t_value = 1.96 + 2.718 / (n - 2); // interpolation for t-student quantile
    // const double t_value = t_student_quantiles[n-11];
    return { m, 2 * t_value * standard_err };
};