#include "CrossSections.hpp"
#include "MonteCarlo.hpp"
#include "MolMass.hpp"
#include "MeanData.hpp"
#include "EnergyData.hpp"
#include "FluxData.hpp"
#include "BulkData.hpp"
#include "ReactionRates.hpp"

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
#include <random>
#include <algorithm>

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
    // E_in_eV in eV for one electron.
    double abs_v = std::sqrt(std::inner_product(v.cbegin(), v.cend(), v.cbegin(), 0.0));
    double E_in_eV =  0.5 * me *  abs_v * abs_v / q0;
    return std::make_pair(abs_v,E_in_eV);
}

void MonteCarlo::maximalCollFreq(){
    // calculates maximal collision rate for gas mixture (in s^-1)
    nu_max = 0.0;

    // Extract total_cross section data from Xsec object:
    const std::vector<double> & energy = Xsec.get_energy();
    const std::vector<double> sigma_tot = Xsec.compute_total_Xsection(mix);

    // Check if the size is coherent:
    if(energy.size() != sigma_tot.size()){
        std::cerr << "Error: energy and total_xsection vectors have different sizes!" << std::endl;
        return;
    }

    // Calculate the maximal collision frequency:
    for(size_t i = 0; i < energy.size(); i++){
        double nu = N * sigma_tot[i] * std::sqrt(2.0 * energy[i] * q0 / me); // s^-1
        if(nu > nu_max) nu_max = nu;
    }
}

void MonteCarlo::initialParticles(){

    // Initialize mean values
    t = {0.0};
    MeanData m(pos_xyz,sigma_xyz,N0);
    
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

    // Update mean vector for the new time step
    mean.emplace_back(mean.back().get_particles(), r[ELECTRONS], v, ne);

    /* If energy update is better here than in "mean" constructor
    energy = std::accumulate(v.cbegin(), v.cend(), 0.0, 
                            [](double sum, const std::array<double, 3>& vi) { return sum + velocity2energy(vi).second; }) / ne;
    mean.back().set_energy(energy);*/
}

void MonteCarlo::update_count_sst(){
    if(T_sst > 0){
        size_t i = 0;
        count_sst = 0;
        for( const double & time : t){
            i++;
            if(time >= T_sst){ count_sst++ };
        }
    }
}

// VERIFY IF THIS VERSION IS BETTER (correct iff every for loop in main only introduces one time step - not sure about that)
// void MonteCarlo::update_count_sst(){
//     if(T_sst > 0){
//          count_sst++;
//     }
// }

void MonteCarlo::updateEnergyData(){
    if( count_sst > 10 ){
        size_t ne = v.size();

        t_total += dt * ne; // sum of all times for all electrons
        std::vector<double> E_in_eV;
        E_in_eV.reserve(ne);
        std::transform(v.begin(), v.end(), std::back_inserter(E_in_eV), [this](const std::array<double, 3>& vi) {
            return velocity2energy(vi).second;
        });

        E.update_energy(E_in_eV,dt,ne,t_total);
        E.compute_distribution_function();
    }    
}

void MonteCarlo::updateFluxData(){
    if(count_sst > 10){
        flux.compute_drift_velocity(v_int,t_total);
        flux.compute_diffusion_const(r[ELECTRONS],v,N);
    }
}

void MonteCarlo::updateBulkData(){
    if (count_sst > 10){
        bulk.update_bulk(t,count_sst,mean,N);
    }
}

void MonteCarlo::updateReactionRates(){

    if (count_sst > 10){

        // 1. REACTION RATES BY COUNTING:
        rates_count.setTime(t,count_sst);
        rates_count.setParticles(mean,count_sst);
        rates_count.computeRates();

        // 2. REACTION RATES BY CONVOLUTION:
        rates_conv.setEnergy(E);
        rates_conv.computeRates();
    
    }   
};

void MonteCarlo::updateCollisionMatrix(){
    // Extract number of electrons:
    int num_particles = v.size();

    // From electron velocity compute energies in eV and extract the velocity modules:
    std::vector<double> E_in_eV(num_particles);
    std::vector<double> v_abs(num_particles);
    std::transform(v.begin(), v.end(), std::back_inserter(E_in_eV), [this](const std::array<double, 3>& vi) {
            return velocity2energy(vi).second;
        });
    std::transform(v.begin(), v.end(), std::back_inserter(v_abs), [this](const std::array<double, 3>& vi) {
            return velocity2energy(vi).first;
        });
    
    // Build collision matrix and compute indeces:
    C.collisionmatrix( num_particles, E_in_eV, v_abs, Xsec.get_energy(), XS.get_full_xs_data(), mix, mgas, E_max, nu_max, N);
    // Update total number of collisions:
    collisions += C.get_collisions();
}

void MonteCarlo::performCollision(const std::string & type){
    // Perform collisions according to the type of collision
    const std::vector<size_t> & ind = C.get_ind(type);
    if(ind.empty()) return; // in case no collisions occoured this time step

    const std::vector<double> & Mass = C.get_mass();
    const std::vector<double> & Loss = C.get_loss();

    if(type == "ELASTIC"){
        elasticCollision(ind, Mass);
    }
    else if(type == "EXCITATION"){
        inelasticCollision(ind, Loss);
    }
    else if(type == "IONIZATION"){
        ionizationCollision(ind, Loss);
    }
    else if(type == "ATTACHMENT"){
        attachmentCollision(ind);
    }
}

std::array<double, 3> MonteCarlo::cross_product(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    // Function to compute the cross product of two 3D vectors
    return {
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
}

void MonteCarlo::elasticCollision(const std::vector<size_t> & ind, const std::vector<double> & Mass){

    double E_1 = 0.0;                       // total energy before collision
    double E_2 = 0.0;                       // total energy after collision
    std::array<double,3> e_x = {1.0, 0.0, 0.0};  // x-direction versor

    std::vector<double> sin_phi(ind.size(),0.0);
    std::vector<double> cos_phi(ind.size(),0.0);
    std::vector<double> sin_xsi(ind.size(),0.0);
    std::vector<double> cos_xsi(ind.size(),0.0);
    std::vector<double> sin_theta(ind.size(),0.0);
    std::vector<double> cos_theta(ind.size(),0.0);

    for(size_t i = 0, i < ind.size(); i++){

        size_t el_index = ind[i];       // electron index
        std::pair<double,double> v2e = velocity2energy(v[el_index]);
        // Compute energy before collision
        E_1 += v2e.second;

        // Compute incident direction of scattering electrons:
        std::array<double,3> e_1 = {v[el_index][0]/v2e.first, v[el_index][1]/v2e.first, v[el_index][2]/v2e.first};

        // Randomly generate phi: azimuthal angle
        double phi = 2 * M_PI * random();
        sin_phi[i] = std::sin(phi);
        cos_phi[i] = std::cos(phi);

        // Randomly generate xsi: electron scattering angle 
        if(isotropic) cos_xsi[i] = 1 - 2 * random();
        else cos_xsi[i] = (2 + v2e.second - 2 * (1+v2e.second) * random()) / v2e.second;
        sin_xsi[i] = std::sqrt(1 - cos_xsi[i] * cos_xsi[i]);

        // Compute theta: angle between x-axis and incident velocity:
        cos_theta[i] = e_1[0];
        sin_theta[i] = std::sqrt(1 - cos_theta[i] * cos_theta[i]);

        // Compute the new direction e_2 of the scattered electron:
        std::array<double,3> cross1 = cross_product(e_1, e_x);
        std::array<double,3> cross2 = cross_product(e_x, e_1);
        std::array<double,3> cross3 = cross_product(e_1, cross2);
        std::array<double, 3> e2;

        // ( Avoid division by zero in case theta is very small):
        double inverse_sin_theta = (sin_theta[i] > 1e-10) ? 1.0 / sin_theta[i] : 0.0;
        for (int j = 0; j < 3; j++) {
            e_2[j] = cos_xsi[i] * e_1[j] +
                    sin_xsi[i] * sin_phi[i] * inverse_sin_theta * cross1[j] +
                    sin_xsi[i] * cos_phi[i] * inverse_sin_theta * cross3[j];
        }

        // Normalize e_2:
        double norm = std::sqrt(e_2[0]*e_2[0] + e_2[1]*e_2[1] + e_2[2]*e_2[2]);
        for (int j = 0; j < 3; ++j) e_2[j] /= norm;

        // Compute energy after the elastic collision:
        E_2 += std::max(0.0, v2e.second*(1 - 2*me/Mass[el_index] * (1 - cos_xsi[i])));

        // Update velocity:
        double v2_abs = std::sqrt(2.0 * E_2 * q0 / me);
        for (int j = 0; j < 3; j++) {
            v[el_index][j] = v2_abs * e_2[j];
        }
    }

    // Update total energy loss counter:
    EnergyLossElastic += E_1 - E_2;
}

void MonteCarlo::inelasticCollision(const std::vector<size_t> & ind, const std::vector<double> & Loss){

    double E_1 = 0.0;                       // total energy before collision
    double E_2 = 0.0;                       // total energy after collision
    std::array<double,3> e_x = {1.0, 0.0, 0.0};  // x-direction versor

    std::vector<double> sin_phi(ind.size(),0.0);
    std::vector<double> cos_phi(ind.size(),0.0);
    std::vector<double> sin_xsi(ind.size(),0.0);
    std::vector<double> cos_xsi(ind.size(),0.0);
    std::vector<double> sin_theta(ind.size(),0.0);
    std::vector<double> cos_theta(ind.size(),0.0);

    for(size_t i = 0, i < ind.size(); i++){

        size_t el_index = ind[i];       // electron index
        std::pair<double,double> v2e = velocity2energy(v[el_index]);
        // Compute energy before collision
        E_1 += v2e.second;

        // Compute incident direction of scattering electrons:
        std::array<double,3> e_1 = {v[el_index][0]/v2e.first, v[el_index][1]/v2e.first, v[el_index][2]/v2e.first};

        // Randomly generate phi: azimuthal angle
        double phi = 2 * M_PI * random();
        sin_phi[i] = std::sin(phi);
        cos_phi[i] = std::cos(phi);

        // Randomly generate xsi: electron scattering angle 
        if(isotropic) cos_xsi[i] = 1 - 2 * random();
        else cos_xsi[i] = (2 + v2e.second - 2 * (1+v2e.second) * random()) / v2e.second;
        sin_xsi[i] = std::sqrt(1 - cos_xsi[i] * cos_xsi[i]);

        // Compute theta: angle between x-axis and incident velocity:
        cos_theta[i] = e_1[0];
        sin_theta[i] = std::sqrt(1 - cos_theta[i] * cos_theta[i]);

        // Compute the new direction e_2 of the scattered electron:
        std::array<double,3> cross1 = cross_product(e_1, e_x);
        std::array<double,3> cross2 = cross_product(e_x, e_1);
        std::array<double,3> cross3 = cross_product(e_1, cross2);
        std::array<double, 3> e2;

        // ( Avoid division by zero in case theta is very small):
        double inverse_sin_theta = (sin_theta[i] > 1e-10) ? 1.0 / sin_theta[i] : 0.0;
        for (int j = 0; j < 3; j++) {
            e_2[j] = cos_xsi[i] * e_1[j] +
                sin_xsi[i] * sin_phi[i] * inverse_sin_theta * cross1[j] +
                sin_xsi[i] * cos_phi[i] * inverse_sin_theta * cross3[j];
        }

        // Normalize e_2:
        double norm = std::sqrt(e_2[0]*e_2[0] + e_2[1]*e_2[1] + e_2[2]*e_2[2]);
        for (int j = 0; j < 3; ++j) e_2[j] /= norm;

        // Compute energy after the elastic collision:
        E_2 += std::max(0.0, v2e.second - Loss[el_index]);

        // Update velocity:
        double v2_abs = std::sqrt(2.0 * E_2 * q0 / me);
        for (int j = 0; j < 3; j++) {
            v[el_index][j] = v2_abs * e_2[j];
        }
    }    

    // Update total energy loss counter:
    EnergyLossInelastic += E_1 - E_2;
}

void MonteCarlo::ionizationCollision(const std::vector<size_t> & ind, const std::vector<double> & Loss){

    double E_1 = 0.0;                       // total energy before collision
    double E_2 = 0.0;                       // total energy after collision
    std::array<double,3> e_x = {1.0, 0.0, 0.0};  // x-direction versor

    // Number of created electrons by ionization:
    int delta_Ne = ind.size();

    std::vector<double> sin_phi(ind.size(),0.0);
    std::vector<double> cos_phi(ind.size(),0.0);
    std::vector<double> sin_xsi(ind.size(),0.0);
    std::vector<double> cos_xsi(ind.size(),0.0);
    std::vector<double> sin_theta(ind.size(),0.0);
    std::vector<double> cos_theta(ind.size(),0.0);

    for(int i = 0, i < delta_Ne; i++){

        size_t el_index = ind[i];       // electron index
        std::pair<double,double> v2e = velocity2energy(v[el_index]);
        // Compute energy before collision
        E_1 += v2e.second;

        // Compute incident direction of scattering electrons:
        std::array<double,3> e_1 = {v[el_index][0]/v2e.first, v[el_index][1]/v2e.first, v[el_index][2]/v2e.first};

        // Randomly generate phi: azimuthal angle
        double phi = 2 * M_PI * random();
        sin_phi[i] = std::sin(phi);
        cos_phi[i] = std::cos(phi);

        // Randomly generate xsi: electron scattering angle 
        if(isotropic) cos_xsi[i] = 1 - 2 * random();
        else cos_xsi[i] = (2 + v2e.second - 2 * (1+v2e.second) * random()) / v2e.second;
        sin_xsi[i] = std::sqrt(1 - cos_xsi[i] * cos_xsi[i]);

        // Compute theta: angle between x-axis and incident velocity:
        cos_theta[i] = e_1[0];
        sin_theta[i] = std::sqrt(1 - cos_theta[i] * cos_theta[i]);

        // Compute the new direction e_2 of the scattered electron:
        std::array<double,3> cross1 = cross_product(e_1, e_x);
        std::array<double,3> cross2 = cross_product(e_x, e_1);
        std::array<double,3> cross3 = cross_product(e_1, cross2);
        std::array<double, 3> e2;

        // ( Avoid division by zero in case theta is very small):
        double inverse_sin_theta = (sin_theta[i] > 1e-10) ? 1.0 / sin_theta[i] : 0.0;
        for (int j = 0; j < 3; j++) {
            e_2[j] = cos_xsi[i] * e_1[j] +
                    sin_xsi[i] * sin_phi[i] * inverse_sin_theta * cross1[j] +
                    sin_xsi[i] * cos_phi[i] * inverse_sin_theta * cross3[j];
        }

        // Normalize e_2:
        double norm = std::sqrt(e_2[0]*e_2[0] + e_2[1]*e_2[1] + e_2[2]*e_2[2]);
        for (int j = 0; j < 3; ++j) e_2[j] /= norm;

        // Compute energy after the elastic collision:
        E_2 += std::max(0.0, v2e.second - Loss[el_index]);

        // Update velocity:
        double v2_abs = std::sqrt(2.0 * W * E_2 * q0 / me); // W is the energy sharing in ionizing collision
        for (int j = 0; j < 3; j++) {
            v[el_index][j] = v2_abs * e_2[j];
        }

        // Add the new electron to the simulation:
        double v_abs_newe = std::sqrt(2.0 * (1-W) * E_2 * q0 / me);
        v.emplace_back(- v_abs_newe * e_2[0], - v_abs_newe * e_2[1], - v_abs_newe * e_2[2]);
        r[ELECTRONS].emplace_back(r[ELECTRONS][el_index][0], r[ELECTRONS][el_index][1], r[ELECTRONS][el_index][2]);

        // Add the new cation:
        r[CATIONS].emplace_back(r[ELECTRONS][el_index][0], r[ELECTRONS][el_index][1], r[ELECTRONS][el_index][2]);
    }

    // If required, enforce electron population conservation:
    if(conserve){
        for( int i = 0; i < delta_Ne; i++){
            // select a random electron:
            int random_index = static_cast<int>(random() * (r[ELECTRONS].size() - 1));
            // remove it from the simulation:
            r[ELECTRONS].erase(r[ELECTRONS].begin() + random_index);
            v.erase(v.begin() + random_index);
        }
    }

    // Add new electrons and cations to mean data
    mean.back().add_new_particles({delta_Ne, delta_Ne, 0});
    // Update total energy loss counter:
    EnergyLossIonization += E_1 - E_2;
}

void MonteCarlo::attachmentCollision(const std::vector<size_t> & ind){

    // Number of created anions / removed electrons by attachment:
    int delta_Ne = ind.size();

    // Note: "ind" is built to be sorted in ascending order (see CollsionData::find_collision_indeces).
    // To remove safely the electrons without changing the other indeces, the loop is performed in reverse order.
    for (auto it = ind.rbegin(); it != ind.rend(); it++) {
    
        size_t el_index = *it;       // electron index

        // Add the new anion to the simulation:
        r[ANIONS].emplace_back(r[ELECTRONS][el_index][0], r[ELECTRONS][el_index][1], r[ELECTRONS][el_index][2]);

        // Remove the electron from the simulation:
        r[ELECTRONS].erase(r[ELECTRONS].begin() + el_index);
        v.erase(v.begin() + el_index);

        // If required, enforce electron population conservation:
        if(conserve){
            // select a random electron:
            int random_index = static_cast<int>(random() * (r[ELECTRONS].size() - 1));
            // clone it to compensate the removed one:
            r[ELECTRONS].push_back(r[ELECTRONS][random_index]);
            v.push_back(v[random_index]);
        }
    }
    
    // Remove electrons and add anions to mean data
    mean.back().add_new_particles({- delta_Ne, 0, delta_Ne});
}

void MonteCarlo::checkSteadyState(){
    if(T_sst == 0.0 && collisions/1e6 >= line && collisions >= col_equ){ 
       
        // check if the interval 80-90 % of energy data is larger than the interval 90-100%:
        int N = mean.size();
        size_t n = std::round( N / 10.0);

        double sum1 = 0.0, sum2 = 0.0;
        for (size_t i = N - 2 * n; i < N - n; ++i) {
            sum1 += mean[i].energy;
        }
        for (size_t i = N - n; i < N; ++i) {
            sum2 += mean[i].energy;
        }

        if(sum1 >= sum2){
            T_sst = t.back();       // last recorded time
            counter = 0;
            collisions = 0;
            line = 1;
        }
    }
}

void MonteCarlo::endSimulation() {
    // End simulation if too many electrons
    if (!conserve && v.size() > Ne_max) {
        converge = 1;
        End = true;
        std::cout << "Simulation ended: maximum number of electrons reached\n";
        return;
    }

    // End simulation if no electrons
    if (!conserve && v.empty()) {
        converge = 2;
        End = true;
        std::cout << "Simulation ended: number of electrons is zero\n";
        return;
    }

    // End simulation if relative errors in w and D are below thresholds
    if (!(bulk.is_empty())) {
        const std::array<double, 3> & w = bulk.get_w();
        const std::array<double, 3> & w_err = bulk.get_w_err();
        const std::array<double, 3> & DN = bulk.get_DN();
        const std::array<double, 3> & DN_err = bulk.get_DN_err();  
        if (std::abs(w_err[2] / w[2]) < std::abs(w_err) &&
            std::abs(DN_err[2] / DN[2]) < std::abs(DN_err)) {

            converge = 0;
            End = true;
            std::cout << "Simulation ended: errors in w < " 
                      << 100 * w_err << "% and D < "
                      << 100 * DN_err << "%\n";
            return;
        }
    }

    // End simulation if number of collisions exceeds maximum
    if (collisions > col_max) {
        converge = 3;
        End = true;
        std::cout << "Simulation ended: maximum number of collisions reached\n";
        return;
    }
}

void MonteCarlo::printOnScreen() {
    if ((collisions / 1e6) >= line) {
        line += 1;

        if (bulk.has_value()) {
            const auto& b = bulk.value();
            const auto& f = flux.value();
            const auto& r = rates;

            std::printf(
                " Werr: %i"
                " DNerr %i"
                " collisions: %i"
                " electrons: %i"
                " E: %.3e eV"
                " w_bulk: %.3e m/s"
                " w_flux: %.3e m/s"
                " DN_bulk: %.2e (ms)^-1"
                " DN_flux: %.2e (ms)^-1"
                " Reff_count: %.2e m^3/s"
                " Reff_calc: %.2e m^3/s"
                " Alpha: %.3e m^-1"
                " Eta: %.3e m^-1\n",

                static_cast<int>(std::abs(b.w_err[2] / b.w[2])),
                static_cast<int>(std::abs(b.DN_err[2] / b.DN[2])),
                collisions,
                mean.back().particles[0], // electrons
                E.E_mean,
                b.w[2],
                f.w[2],
                b.DN[2],
                f.DN[2],
                r.count.eff,
                r.conv.eff,
                r.count.ion_tot * 2.4e25 / b.w[2],
                r.count.att_tot * 2.4e25 / b.w[2]
            );

        } else {
            std::printf(
                " collisions: %i"
                " electrons: %i"
                " mean energy: %.2e\n",
                collisions,
                mean.back().particles[0],
                mean.back().energy
            );
        }
    }
}
