#ifndef MEANDATA_H
#define MEANDATA_H

#include <array>
#include <vector>
#include <iostream>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "Common.hpp"

// Object MeanData stores the mean data of the system at a given time
class MeanData
{
    public:

    // Default constructor:
    MeanData() = default;
    // Constructor for initializing electrons:
    MeanData(const std::array<double,3> & r, const std::array<double,3> & s, const int ne)
        : particles({ne, 0, 0}),
          position(r),
          var({s[0]*s[0], s[1]*s[1], s[2]*s[2]}), // variance is sigma^2
          velocity({0.0, 0.0, 0.0}),
          energy(0.0) {};
    // Constructor that computes mean data from input matrices:
    MeanData(const std::array<int,mc::PARTICLES_TYPES> & p, const mc::MATRIX & r, const mc::MATRIX & v)
        : particles(p){

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


    // Add particles (e.g. created by ionization):
    void add_new_particles(const std::array<int,mc::PARTICLES_TYPES> & p){
        particles[mc::ELECTRONS] += p[mc::ELECTRONS];
        particles[mc::CATIONS] += p[mc::CATIONS];
        particles[mc::ANIONS] += p[mc::ANIONS];
    };

    // Setters:
    void set_particles(const std::array<int,mc::PARTICLES_TYPES> & p){
        particles = p;
    };
    void set_position(const std::array<double,3> & pos){
        position = pos;
    };
    void set_variance(const std::array<double,3> & var){
        this->var = var;
    };
    void set_velocity(const std::array<double,3> & vel){
        velocity = vel;
    };
    void set_energy(const double & e){
        energy = e;
    }

    // Getters:
    const std::array<int,mc::PARTICLES_TYPES> & get_particles() const{
        return particles;
    }
    int get_N_electrons() const{
        return particles[mc::ELECTRONS];
    }
   const  std::array<double,3> & get_position() const{
        return position;
    }
    const std::array<double,3> & get_variance() const{
        return var;
    }
    const std::array<double,3> & get_velocity() const{
        return velocity;
    }
    double get_energy() const{
        return energy;
    }

    private:

    std::array<int,mc::PARTICLES_TYPES> particles; // # of particles per each type
    
    // Electrons mean data:
    std::array<double, 3> position;
    std::array<double, 3> var;
    std::array<double, 3> velocity;
    double energy;
};

#endif