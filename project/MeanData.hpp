#ifndef MEANDATA_H
#define MEANDATA_H

#include <array>
#include <vector>
#include <iostream>
#include <numeric>
#include <cmath>
#include <algorithm>

enum ParticleType {ELECTRONS = 0, CATIONS, ANIONS, PARTICLES_TYPES};
typedef std::vector<std::array<double,3>> MATRIX;

// Object MeanData stores the mean data of the system at a given time
class MeanData
{
    public:

    // Default constructor:
    MeanData() = default;
    // Constructor for initializing electrons:
    MeanData(const std::array<double,3> & r, const std::array<double,3> & s, const int & ne)
        : particles({ne, 0, 0}),
          position(r),
          sigma(s),
          energy(0.0),
          velocity({0.0, 0.0, 0.0}) {};
    // Constructor that computes mean data from input matrices:
    MeanData(const std::array<int,PARTICLES_TYPES> & p, const MATRIX & r, const MATRIX & v, const double & ne)
        : particles(p){
        
            // Control that number of electrons is non-zero.
            if (ne == 0) {
                std::cerr << "Warning: ne is zero. Mean values will be zero.\n";
                return;
            }
        
            // Compute mean position:
            for (size_t i = 0; i < 3; i++) {
                position[i] = std::accumulate(r.cbegin(), r.cend(), 0.0, 
                                [i](double sum, const std::array<double, 3>& ri) { return sum + ri[i]; }) / ne;
            }
            // Compute mean velocity:
            for (size_t i = 0; i < 3; i++) {
                velocity[i] = std::accumulate(v.cbegin(), v.cend(), 0.0, 
                                [i](double sum, const std::array<double, 3>& vi) { return sum + vi[i]; }) / ne;
            }

            // Compute mean energy:
            double me = 9.10938291e-31; // electron mass
            double q0 = 1.60217657e-19; // electron charge
            energy = std::accumulate(v.cbegin(), v.cend(), 0.0, 
                           [me, q0](double sum, const std::array<double, 3>& vi){
                                double abs_v = std::sqrt(std::inner_product(vi.cbegin(), vi.cend(), vi.cbegin(), 0.0));
                                return sum + 0.5 * me * abs_v * abs_v / q0;
                           }) / ne;

            // Compute standard deviation:
            for (size_t i = 0; i < 3; i++) {
                sigma[i] = std::sqrt(std::accumulate(r.cbegin(), r.cend(), 0.0,
                                [i, this](double sum, const std::array<double, 3>& ri) { 
                                    return sum + std::pow(ri[i] - position[i], 2);
                                }) / ne);
            }
    };

    // Add particles (e.g. created by ionization):
    void add_new_particles(std::array<int,PARTICLES_TYPES> & p){
        particles[ELECTRONS] += p[ELECTRONS];
        particles[CATIONS] += p[CATIONS];
        particles[ANIONS] += p[ANIONS];
    };

    // Setters:
    void set_particles(const std::array<int,PARTICLES_TYPES> & p){
        particles = p;
    };
    void set_position(const std::array<double,3> & pos){
        position = pos;
    };
    void set_sigma(const std::array<double,3> & sig){
        sigma = sig;
    };
    void set_velocity(const std::array<double,3> & vel){
        velocity = vel;
    };
    void set_energy(const double & e){
        energy = e;
    }

    // Getters:
    const std::array<int,PARTICLES_TYPES> & get_particles() const{
        return particles;
    }
    int get_N_electrons() const{
        return particles[ELECTRONS];
    }
   const  std::array<double,3> & get_position() const{
        return position;
    }
    const std::array<double,3> & get_sigma() const{
        return sigma;
    }
    const std::array<double,3> & get_velocity() const{
        return velocity;
    }
    double get_energy() const{
        return energy;
    }

    private:

    std::array<int,PARTICLES_TYPES> particles; // # of particles per each type
    
    // Electrons mean data:
    std::array<double, 3> position;
    std::array<double, 3> sigma;
    std::array<double, 3> velocity;
    double energy;
};

#endif