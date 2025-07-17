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
    void compute_drift();
    void compute_diffusion(const double N);
    const std::array<double,2> linear_regression(const std::vector<double>& y) const;

    //----------------------------------------------------------------------------------------------------------//
    //-------------------------------------- FOR DEBUGGING PURPOSES --------------------------------------------//
    const std::array<double,143> t_student_quantiles = {
        1.859548, 1.833113, 1.812461, 1.795885, 1.782288, 1.770933, 1.76131, 1.75305,
        1.745884, 1.739607, 1.734064, 1.729133, 1.724718, 1.720743, 1.717144, 1.713872,
        1.710882, 1.708141, 1.705618, 1.703288, 1.701131, 1.699127, 1.697261, 1.695519, 1.693889,
        1.69236, 1.690924, 1.689572, 1.688298, 1.687094, 1.685954, 1.684875, 1.683851, 1.682878, 1.681952,
        1.681071, 1.68023, 1.679427, 1.67866, 1.677927, 1.677224, 1.676551, 1.675905, 1.675285, 1.674689, 1.674116,
        1.673565, 1.673034, 1.672522, 1.672029, 1.671553, 1.671093, 1.670649, 1.670219, 1.669804, 1.669402, 1.669013,
        1.668636, 1.668271, 1.667916, 1.667572, 1.667239, 1.666914, 1.6666, 1.666294, 1.665996, 1.665707, 1.665425,
        1.665151, 1.664885, 1.664625, 1.664371, 1.664125, 1.663884, 1.663649, 1.66342, 1.663197, 1.662978, 1.662765,
        1.662557, 1.662354, 1.662155, 1.661961, 1.661771, 1.661585, 1.661404, 1.661226, 1.661052, 1.660881, 1.660715,
        1.660551, 1.660391, 1.660234, 1.660081, 1.65993, 1.659782, 1.659637, 1.659495, 1.659356, 1.659219, 1.659085,
        1.658953, 1.658824, 1.658697, 1.658573, 1.65845, 1.65833, 1.658212, 1.658096, 1.657982, 1.65787, 1.657759,
        1.657651, 1.657544, 1.657439, 1.657336, 1.657235, 1.657135, 1.657037, 1.65694, 1.656845, 1.656752, 1.656659,
        1.656569, 1.656479, 1.656391, 1.656305, 1.656219, 1.656135, 1.656052, 1.65597, 1.65589, 1.655811, 1.655732,
        1.655655, 1.655579, 1.655504, 1.65543, 1.655357, 1.655285, 1.655215, 1.655145, 1.655076
    };
    //----------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------------------------------------------------//
};

#endif