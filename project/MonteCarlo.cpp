#include "cross_s.hpp"
#include "MonteCarlo.hpp"
#include "MolMass.hpp"
#include "Mesh.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <cmath>
#include <utility>

void MonteCarlo::checkFractionSum(){
    double sumMix = std::accumulate(this->mix.cbegin(),this->mix.cend(),0.0);
    
    // checks if sum of gas fractions is equal to 1.
    // If not the case: last entry of mix will be corrected.               

    if(sumMix != 1){
            this->mix.back() = 1 - std::accumulate(this->mix.cbegin(),this->mix.cend()-1,0.0);
            std::cout << "Sum of mixing ratios (in mix) NOT equal to one!" << std::endl;
    }
}

void MonteCarlo::mass_in_kg(){
    // Computes mass of the gas species in kg
    this->mgas.clear(); // In case the method is called more than once by mistake
    MolMass mol_mass;
    for(std::string & sub : this->gas){
        mol_mass.set_substance(sub);
        mol_mass.Compute_M();
        for( const double & M : mol_mass.get_M()){
            this->mgas.emplace_back(M/(this->Na * 1000));
        }
    }
}

void MonteCarlo::gasNumberDensity(){
    // Calculates the gas number density N (in m^-3) by the ideal gas law
    // from pressure p and temperature Temp.
    this->N = this->p / (this->kB * this->Temp);
}

std::pair<double,double> MonteCarlo::velocity2energy_in_ev(const std::vector<double> & v){
    // calculates absolute value of velocity abs_v and energy
    // E_in_eV in eV.
    double abs_v = std::sqrt(std::inner_product(v.cbegin(), v.cend(), v.cbegin(), 0.0));
    double E_in_eV =  0.5 * this->me *  abs_v * abs_v / this->q0;
    return std::make_pair(abs_v,E_in_eV);
}

// void MonteCarlo::maximalCollFreq(){}

void MonteCarlo::initialParticles(){

    // Initialize mean values
    this->t = {0.0};
    this->mean.energy = 0.0;
    this->mean.position = this->pos_xyz; // pos_xyz defined by user
    this->mean.sigma = this->sigma_xyz;  // sigma_xyz defined by user
    this->mean.velocity = {0.0, 0.0, 0.0};

    // Initialize number of particles
    this->mean.particles[ELECTRONS] = this->N0;
    this->mean.particles[CATIONS] = 0;
    this->mean.particles[ANIONS] = 0;

    // Allocate memory to position, velocity and acceleration vectors for electrons
    this->r[ELECTRONS].resize(this->N0);
    this->v.resize(this->N0);
    this->a.resize(this->N0);

    // Random number generators
    std::default_random_engine gen;
    std::normal_distribution<double> randn(0.0,1.0); // standard Gaussian distribution
    std::uniform_real_distribution<> randu(0.0,1.0); // [0,1] uniform distribution
    // Should not be better a randu(-1.0,1.0) ?
    
    // Initialize electrons position
    for( auto &re : this->r[ELECTRONS]){
        for( size_t i=0; i<=2; i++){ 
            // Initialize to pos_xyz + noise if sigma_xyz != 0
            re[i] = this->pos_xyz[i] + this->sigma_xyz[i] * randn(gen); 
        } 
    }
    
    // Initialize electrons velocities
    for( auto &vel : this->v){
        for( size_t i=0; i<=2; i++){
            
            vel[i] = 0.0; // Initialize electrons with zero velocity

            // Add noise to avoid singularities in MonteCarlo::elasticCollision 
            vel[i] += 1e-6 * randu(gen);

            // MATLAB CODE NOT CLEAR HERE:
            // Apparently it implements noise, in fact it is not added because of the
            // coefficient "max(max(abs(v)))" which is ALWAYS 0 (initialized literally one line above!).
            // Here I ignored such coefficient but I am not sure this is correct.
        } 
    }
}

void MonteCarlo::surfaceInteraction(){
    auto it_r = this->r[ELECTRONS].begin(); 
    auto it_v = this->v.begin();
    // loop over electrons
    while(it_r != this->r[ELECTRONS].end()){
        auto &re = *it_r;
        if( re[0] < 0 || re[0] > this->Lx ||
            re[1] < 0 || re[1] > this->Ly ||
            re[2] < 0 || re[2] > this->Lz ){
                // remove electrons outside the boundary
                it_r = this->r[ELECTRONS].erase(it_r);
                it_v = this->v.erase(it_v);
        } 
        else{
            it_r++;
            it_v++;
        }
    }
}

//void MonteCarlo::particle2density(){
    // calculates electron number density rho (in 1/m^3)

    // Implement these definition + interpolation method in class Mesh
    //double dx = this->x[1] - this->x[0];
    //double dy = this->y[1] - this->y[0];
    //double dz = this->z[1] - this->z[0];

    // Allocate memory for density vector
    //this->rho.reserve(this->r[ELECTRONS].size);

    //for( auto &rr : this->r ){
        // x = mesh.interpolation(rr[1]);
        // y = ...
        // z = ...
    //}        

    // 
//}


