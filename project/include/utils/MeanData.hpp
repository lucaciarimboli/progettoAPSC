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

    // Constructors:
    MeanData(const std::array<double,3> & r, const std::array<double,3> & s, const int ne):
        particles({ne, 0, 0}),
        position(r),
        var({s[0]*s[0], s[1]*s[1], s[2]*s[2]}), // variance is sigma^2
        velocity({0.0, 0.0, 0.0}),
        energy(0.0)
    {};
    MeanData(const std::array<int,mc::PARTICLES_TYPES> & p, const mc::MATRIX & r, const mc::MATRIX & v);

    // Public Method:
    void add_new_particles(const std::array<int,mc::PARTICLES_TYPES> & p);

    // Setters:
    void set_particles(const std::array<int,mc::PARTICLES_TYPES> & p){ particles = p; };
    void set_position(const std::array<double,3> & pos){ position = pos; };
    void set_variance(const std::array<double,3> & var){ this->var = var; };
    void set_velocity(const std::array<double,3> & vel){ velocity = vel; };
    void set_energy(const double & e){ energy = e; };

    // Getters:
    const std::array<int,mc::PARTICLES_TYPES> & get_particles() const{ return particles; }
    int get_N_electrons() const{ return particles[mc::ELECTRONS]; }
    const  std::array<double,3> & get_position() const{ return position; }
    const std::array<double,3> & get_variance() const{ return var; }
    const std::array<double,3> & get_velocity() const{ return velocity; }
    double get_energy() const{ return energy; }

    private:

    // Class Members:
    std::array<int,mc::PARTICLES_TYPES> particles; // # of particles per each type
    std::array<double, 3> position;                // mean position
    std::array<double, 3> var;                     // position variance
    std::array<double, 3> velocity;                // mean velocity
    double energy;                                 // mean energy [eV]
};

#endif