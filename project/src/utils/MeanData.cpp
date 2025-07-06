#include "utils/MeanData.hpp"

MeanData::MeanData(const std::array<double,3> & r, const std::array<double,3> & s, const int ne):
    // Constructor to set initial mean value
    particles({ne, 0, 0}),
    position(r),
    var({s[0]*s[0], s[1]*s[1], s[2]*s[2]}), // variance is sigma^2
    velocity({0.0, 0.0, 0.0}),
    energy(0.0)
{};

MeanData::MeanData(const std::array<int,mc::PARTICLES_TYPES> & p, const mc::MATRIX & r, const mc::MATRIX & v)
    : particles(p){

    // Constructor that computes mean data from input matrices:

    const int ne = r.size();

    // Compute mean position:
    position[0] = std::accumulate(r.cbegin(), r.cend(), 0.0, 
                        [](double sum, const std::array<double, 3>& ri) { return sum + ri[0]; }) / ne;
    position[1] = std::accumulate(r.cbegin(), r.cend(), 0.0, 
                        [](double sum, const std::array<double, 3>& ri) { return sum + ri[1]; }) / ne;
    position[2] = std::accumulate(r.cbegin(), r.cend(), 0.0, 
                        [](double sum, const std::array<double, 3>& ri) { return sum + ri[2]; }) / ne;

    // Compute mean velocity:
    velocity[0] = std::accumulate(v.cbegin(), v.cend(), 0.0, 
                        [](double sum, const std::array<double, 3>& vi) { return sum + vi[0]; }) / ne;
    velocity[1] = std::accumulate(v.cbegin(), v.cend(), 0.0, 
                        [](double sum, const std::array<double, 3>& vi) { return sum + vi[1]; }) / ne;
    velocity[2] = std::accumulate(v.cbegin(), v.cend(), 0.0, 
                        [](double sum, const std::array<double, 3>& vi) { return sum + vi[2]; }) / ne;

    // Compute mean energy:
    energy = std::accumulate(v.cbegin(), v.cend(), 0.0, 
                    [](double sum, const std::array<double, 3>& vi){
                        double abs_v = std::sqrt(std::inner_product(vi.cbegin(), vi.cend(), vi.cbegin(), 0.0));
                        return sum + 0.5 * mc::me * abs_v * abs_v / mc::q0;
                    }) / ne;

    // Compute variance of position:
    var[0] = std::accumulate(r.cbegin(), r.cend(), 0.0,
                        [this](double sum, const std::array<double, 3>& ri) {
                            return sum + (ri[0] - position[0]) * (ri[0] - position[0]);
                        }) / (ne - 1);
    var[1] = std::accumulate(r.cbegin(), r.cend(), 0.0,
                        [this](double sum, const std::array<double, 3>& ri) {
                            return sum + (ri[1] - position[1]) * (ri[1] - position[1]);
                        }) / (ne - 1);
    var[2] = std::accumulate(r.cbegin(), r.cend(), 0.0,
                        [this](double sum, const std::array<double, 3>& ri) {
                            return sum + (ri[2] - position[2]) * (ri[2] - position[2]);
                        }) / (ne - 1);
};

void MeanData::add_new_particles(const std::array<int,mc::PARTICLES_TYPES> & p){
    // Add particles (e.g. created by ionization):
    particles[mc::ELECTRONS] += p[mc::ELECTRONS];
    particles[mc::CATIONS] += p[mc::CATIONS];
    particles[mc::ANIONS] += p[mc::ANIONS];
};

