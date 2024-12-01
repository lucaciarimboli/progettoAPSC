// #include "cross_s.hpp"
#include "prova.hpp"
#include "MeanData.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <cmath>
#include <utility>

std::pair<double,double> PROVA::velocity2energy(const std::array<double,3> & v){
    // calculates absolute value of velocity abs_v and energy
    // E_in_eV in eV.
    double abs_v = std::sqrt(std::inner_product(v.cbegin(), v.cend(), v.cbegin(), 0.0));
    double E_in_eV =  0.5 * me *  abs_v * abs_v / q0;
    return std::make_pair(abs_v,E_in_eV);
}

void PROVA::initialParticles(){

    // Initialize mean values
    t = {0.0};
    MeanData m(pos_xyz,sigma_xyz);
    m.set_particles(N0,0,0);
    
    mean.clear();
    mean.emplace_back(m);

    // Allocate memory to position, velocity and acceleration vectors for electrons
    r[ELECTRONS].resize(N0);
    v.resize(N0);
    a.resize(N0);

    // Random number generators
    std::default_random_engine gen;
    std::normal_distribution<double> randn(0.0,1.0); // standard Gaussian distribution
    std::uniform_real_distribution<> randu(0.0,1.0); // [0,1] uniform distribution
    // Should not be better a randu(-1.0,1.0) ?
    
    // Initialize electrons position
    for( auto &re : r[ELECTRONS]){
        for( size_t i=0; i<=2; i++){ 
            // Initialize to pos_xyz + noise if sigma_xyz != 0
            re[i] = pos_xyz[i] + sigma_xyz[i] * randn(gen); 
        } 
    }
    
    // Initialize electrons velocities
    for( auto &vel : v){
        for( size_t i=0; i<=2; i++){
            
            vel[i] = 0.0; // Initialize electrons with zero velocity

            // Add noise to avoid singularities in MonteCarlo::elasticCollision 
            // vel[i] += 1e-6 * randu(gen);

            // MATLAB CODE NOT CLEAR HERE:
            // Apparently it implements noise, in fact it is not added because of the
            // coefficient "max(max(abs(v)))" which is ALWAYS 0 (initialized literally one line above!).
            // Here I ignored such coefficient but I am not sure this is correct.
        } 
    }
}