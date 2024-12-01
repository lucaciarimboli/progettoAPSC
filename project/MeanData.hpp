#ifndef MEANDATA_H
#define MEANDATA_H

#include <array>

enum ParticleType {ELECTRONS = 0, CATIONS, ANIONS, PARTICLES_TYPES};

// Object MeanData stores the mean data of the system at a given time
class MeanData
{
    public:

    // Default constructor:
    MeanData() = default;
    // Constructor for initializing electrons:
    MeanData(const std::array<double,3> & p, const std::array<double,3> & s)
        : position{p},
          sigma{s},
          energy(0.0),
          velocity({0.0, 0.0, 0.0}) {};

    // Setters:
    void set_particles(const std::array<int,PARTICLES_TYPES> & p){
        particles = p;
    };
    void set_position(const double & x, const double & y, const double & z){
        position = {x,y,z};
    };
    void set_sigma(const double & sx, const double & sy, const double & sz){
        sigma = {sx,sy,sz};
    };
    void set_velocity(const double & vx, const double & vy, const double & vz){
        velocity = {vx,vy,vz};
    };
    void set_energy(const double & e){
        energy = e;
    }

    // Getters:
    std::array<int,PARTICLES_TYPES> & get_particles() const{
        return particles;
    }
    int get_N_electrons() const{
        return particles[ELECTRONS];
    }
    std::array<double,3> & get_position() const{
        return position;
    }
    std::array<double,3> & get_sigma() const{
        return sigma;
    }
    std::array<double,3> & get_velocity() const{
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