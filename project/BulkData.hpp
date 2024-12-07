#ifndef ENERGYDATA_H
#define ENERGYDATA_H

#include <array>
#include <vector>
#include "MeanData.hpp"

class MeanData;

class BulkData {
    public:
    // Constructor
    BulkData() : w{0.0, 0.0, 0.0}, w_err{0.0, 0.0, 0.0},
                 DN{0.0, 0.0, 0.0}, DN_err{0.0, 0.0, 0.0} {}

    void drift(){
        // Each component of y is a vector with the corrisponding component of the position for each time after sst.
        std::array<std::vector<double>,3> y;
        for( size_t j = 0; j < mean.size(); j++){
            std::array<double,3> pos = mean[j].get_position();
            for( size_t i = 0; i < 3; i++){
                y[i][j].emplace_back(pos[i]);
            }
        }
        
        x_max = max(x);

        for(i = 0; i < 3; i++){
            double y_max = max(y);
            auto B = linear_regression(y[j]); // ??????
            w[i] = B.first[i] * y_max / x_max;
            // w_err[i] = 0.25 * diff(B.second[i]) * y_max / x_max;
        }
        
    };

    void diffusion(const double & N){
        std::array<std::vector<double>,3> y;
        for( size_t i = 0; i < mean.size(); i++){
            std::array<double,3> sig = mean[i].get_sigma();
            for( size_t j = 0; j < 3; j++){
                y[j][i].emplace_back(sig[j]);
            }
        }
        std::array<double,3> y_max = max(y);
        // Da completare - simile a drift();
    }

    void set_t(const std::vector<double> & x){ t = x; }
    void set_mean(const std::vector<MeanData> & y){ mean = y; }

    private:
    std::array<double, 3> w;      // Drift velocity
    std::array<double, 3> w_err; // Error in drift velocity
    std::array<double, 3> DN;    // Diffusion constant
    std::array<double, 3> DN_err; // Error in diffusion constant

    std::vector<double> t;        // time starting from T_sst
    std::vector<MeanData> mean;   // mean data of steady states

    // Compute max value of a std::vector<double>
    std::array<double,3> max(std::vector<double> & y){
        double y_max = 0.0;
        for( double & yy : y){
                if(yy > y_max){
                    y_max = yy;
                }
            }
        return y_max;
    };

    // diff function of matlab

    // Linear regression that should work as matlab function "regress".
    // y first input is a vector 
    std::pair<double,double> linear_regression(const std::vector<double>& y, const std::vector<double>& x);
};

#endif