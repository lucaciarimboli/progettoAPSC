// #include "cross_s.hpp"
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
#include <utility>

void PROVA::mass_in_kg(){
    // converts mass (mgas) for the gas species from a.u. into kg.
    this->mgas.clear(); // In case the method is called more than once by mistake
    MolMass mol_mass;
    for(std::string & sub : this->gas){
        mol_mass.set_substance(sub);
        mol_mass.Compute_M();
        std::cout << mol_mass.get_substance() << ": \n" << std::endl;
        for( const double & M : mol_mass.get_M()){
            std::cout << M << "\n" << std::endl;
            this->mgas.emplace_back(M/(this->Na * 1000));
        }
    }
}

void PROVA::gasNumberDensity(){
    // Calculates the gas number density N (in m^-3) by the ideal gas law
    // from pressure p and temperature Temp.
    this->N = this->p / (this->kB * this->Temp);
}

std::pair<double,double> PROVA::velocity2energy_in_ev(const std::vector<double> & v){
    // calculates absolute value of velocity abs_v and energy
    // E_in_eV in eV.
    double abs_v = std::sqrt(std::inner_product(v.cbegin(), v.cend(), v.cbegin(), 0.0));
    double E_in_eV =  0.5 * this->me *  abs_v * abs_v / this->q0;
    return std::make_pair(abs_v,E_in_eV);
}