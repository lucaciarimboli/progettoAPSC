#include "prova.hpp"
#include "MolMass.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <cmath>

int main(){
    
    PROVA prova;
    prova.gas = {"O2", "[Ar,   Au]", "N2(H2O)CO2"};
    std::cout << "Mass in a.u" << "\n" << std::endl;

    prova.mass_in_kg();
    
    std::cout << "Ora in kg" << "\n"  << std::endl;

    for( const auto &s : prova.mgas){
        std::cout << s << "\n" << std::endl;
    }

    return 0;
}

//enum ParticleType {ELECTRONS = 0, CATIONS, ANIONS, PARTICLES_TYPES};
//typedef std::vector<std::array<double,3>> POSITIONS;

//std::array<POSITIONS,PARTICLES_TYPES> r;

// array of initial mean position of initial gaussian distributed electrons in x,y and z direction
//std::array<double,3> pos_xyz  = {0, 0, 0};
// array of initial broadening of initial gaussian distributed electrons in x,y and z direction
//std::array<double,3> sigma_xyz = {0, 0, 0};