#include "cross_s.hpp"
#include "MonteCarlo.hpp"
#include "MolMass.hpp"
#include "MeanData.hpp"
#include "EnergyData.hpp"
#include "FluxData"
#include "BulkData"

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <cmath>
#include <utility>

void MonteCarlo::checkFractionSum(){
    double sumMix = std::accumulate(mix.cbegin(),mix.cend(),0.0);
    
    // checks if sum of gas fractions is equal to 1.
    // If not the case: last entry of mix will be corrected.           

    if(sumMix != 1){
            mix.back() = 1 - std::accumulate(mix.cbegin(),mix.cend()-1,0.0);
            std::cout << "Sum of mixing ratios (in mix) NOT equal to one!" << std::endl;
    }
}

void MonteCarlo::mass_in_kg(){
    // Computes mass of the gas species in kg
    mgas.clear(); // In case the method is called more than once by mistake
    MolMass mol_mass;
    for(std::string & sub : gas){
        mol_mass.set_substance(sub);
        mol_mass.Compute_M();
        for( const double & M : mol_mass.get_M()){
            mgas.emplace_back(M/(Na * 1000));
        }
    }
}

void MonteCarlo::gasNumberDensity(){
    // Calculates the gas number density N (in m^-3) by the ideal gas law
    // from pressure p and temperature Temp.
    N = p / (kB * Temp);
}

std::pair<double,double> MonteCarlo::velocity2energy(const std::array<double,3> & v){
    // calculates absolute value of velocity abs_v and energy
    // E_in_eV in eV.
    double abs_v = std::sqrt(std::inner_product(v.cbegin(), v.cend(), v.cbegin(), 0.0));
    double E_in_eV =  0.5 * me *  abs_v * abs_v / q0;
    return std::make_pair(abs_v,E_in_eV);
}

// void MonteCarlo::maximalCollFreq(){} // USES CLASS "cross_sect"

void MonteCarlo::initialParticles(){

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

void MonteCarlo::surfaceInteraction(){
    auto it_r = r[ELECTRONS].begin(); 
    auto it_v = v.begin();
    // loop over electrons
    while(it_r != r[ELECTRONS].end()){
        auto &re = *it_r;
        if( re[0] < 0 || re[0] > Lx ||
            re[1] < 0 || re[1] > Ly ||
            re[2] < 0 || re[2] > Lz ){
                // remove electrons outside the boundary
                it_r = r[ELECTRONS].erase(it_r);
                it_v = v.erase(it_v);
        } 
        else{
            it_r++;
            it_v++;
        }
    }
}

void MonteCarlo::set_E(){
    // Works for the case: EN constant, uniform and user defined.
    // "solvePoisson()" is an extension of this method.
    E_x = 0;
    E_y = 0;
    E_z = EN * N * 1e-21; // set E // z direction !
}

double MonteCarlo::random(){
        // Generates a random number from a U[0,1]
        std::default_random_engine gen;
        std::uniform_real_distribution<> randu(0.0,1.0); // [0,1] uniform distribution
        return randu(gen);
}

std::vector<double> MonteCarlo::random(const unsigned int N){
    //generates an N-vector of random numbers from a U[0,1]
    std::vector<double,N> p;
    p.resize(N);

    std::default_random_engine gen;
    std::uniform_real_distribution<> randu(0.0,1.0); // [0,1] uniform distribution
    for( double & pp : p){
        pp = randu(gen);
    }
    return p;
}

void MonteCarlo::freeFlight(){
    // performs non-collissional flight for electrons in electric field
    dt = - std::log(random()) / nu_max; // generates time step

    // Update vector time:
    t.emplace_back(t.back() + dt);
    
    double ne = v.size(); // number of electrons
    MATRIX a; // acceleration matrix
    a[0].resize(ne, q0/me * E_x);
    a[1].resize(ne, q0/me * E_y);
    a[2].resize(ne, q0/me * E_z);

    double dt2 = std::pow(dt,2) // to avoid redundant calculations 

    if(T_sst > 0.0){
        v_int.resize(ne);
        v2_int.resize(ne);
        double dt3 = std::pow(dt,3);
        for(size_t i = 0; i < ne; i++){     // i indicates the electron
            for(size_t j = 0; j < 3; j++){  // j indicates the coordinate (x,y,z)
                // integrated velocity
                v_int[i][j] = v[i][j]*dt + 0.5 * a[i][j] * dt2;
                // integrated velocity squared
                v2_int[i][j] = std::pow(v[i][j],2)*dt + a[i][j]*v[i][j] * dt2 + std::pow(a[i][j],2)/3 * dt3;
            }
        }
    }

    // Update space and velocity:
    for(size_t i = 0; i < ne; i++){
            for(size_t j = 0; j < 3; j++){
                r[ELECTRONS][i][j] += v[i][j]*dt + 0.5 * a[i][j] * dt2;
                v[i][j] += a[i][j] * dt;
        }
    }
}

void MonteCarlo::collectMeanData(){
    MeanData m;
    m.set_particles(mean.back().get_particles());

    // Average position, velocity and energy of electrons:
    sum_x = 0.0;
    sum_y = 0.0;
    sum_z = 0.0;
    sum_vx = 0.0;
    sum_vy = 0.0;
    sum_vz = 0.0;
    sum_e = 0.0;
    int ne = m.get_N_el();
    for(size_t i = 0; i < ne; i++){
        sum_x += r[ELECTRONS][i][0];
        sum_y += r[ELECTRONS][i][1];
        sum_z += r[ELECTRONS][i][2];
        sum_vx += v[i][0];
        sum_vy += v[i][1];
        sum_vz += v[i][2];
        sum_e += velocity2energy(v[i]).second;
    }
    m.set_position(sum_x/ne,sum_y/ne,sum_z/ne);
    m.set_velocity(sum_vx/ne,sum_vy/ne,sum_vz/ne);
    m.set_energy(sum_e/ne);

    // Standard deviation
    double sigma_x = 0.0;
    double sigma_y = 0.0;
    double sigma_z = 0.0;
    const auto& mean_position = m.get_position();
    for(const auto & pos : r[ELECTRONS]){
        sigma_x += std::pow(pos[0]-mean_position[0],2);
        sigma_y += std::pow(pos[1]-mean_position[1],2);
        sigma_z += std::pow(pos[2]-mean_position[2],2);
    }
    m.set_sigma(std::sqrt(sigma_x/ne),std::sqrt(sigma_y/ne),std::sqrt(sigma_z/ne));

    // Update mean vector for the new time step
    mean.emplace_back(m);
}

void MonteCarlo::update_flag_sst(){
    if(T_sst > 0){
        unsigned int count = 0;
        for( const double & time : t){
            if(time >= T_sst) count++;
        }
    }
    if(count > 10){
        flag_sst = 1;
    }
}

void MonteCarlo::updateEnergyData(){
    if( flag_sst ){
        size_t ne = v.size();

        t_total += dt * ne;
        std::vector<double> E_in_eV(ne);
        for (size_t i = 0; i < ne; i++) {
                E_in_eV[i] = velocity2energy(v[i]).second;
        }

        E.update_energy(E_in_eV,dt,ne);
        E.set_E_mean(t_total);
        E.compute_distribution_function();
    }    
}

void MonteCarlo::fluxData(){
    if(flag_sst){
        flux.compute_drift_velocity(v_int,t_total);
        flux.compute_diffusion_const(r[ELECTRONS],v,N);
    }
}
