#include <iostream>
#include <vector>
#include <string>
#include <random>
#include "../include/utils/CrossSectionsData.hpp"
#include "../include/Common.hpp"

int main(){
    const std::vector<std::string> species = {"N2", "N", "O2", "O", "H2O", "H2"};
    const size_t n_species = species.size();

    const std::vector<double> mix(n_species,0.0);   // not needed for plot
    const double rho = 0.0;                 // not needed for plot
    const double E_max = 1e5;

    const CrossSectionsData Xsec(species,E_max,mix,rho);

    // Save energy grid values:
    const std::vector<double> energy = Xsec.get_energy();
    const size_t n = energy.size();
    std::ofstream file("plots/Xsec/energy.csv");
    for(size_t i = 0; i<n; i++){
        file << energy[i] << "\n";
    }
    
    for( size_t spec_idx = 0; spec_idx < n_species; spec_idx++){

        std::vector<table> Xsec_eff = Xsec.get_Xsections(spec_idx,mc::EFFECTIVE);
        std::vector<table> Xsec_att = Xsec.get_Xsections(spec_idx,mc::ATTACHMENT);
        std::vector<table> Xsec_exc = Xsec.get_Xsections(spec_idx,mc::EXCITATION);
        std::vector<table> Xsec_ela = Xsec.get_Xsections(spec_idx,mc::ELASTIC);

        // Build vectors with cross sections data for every energy level:
        std::vector<double> eff(n);
        std::vector<double> att(n);
        std::vector<double> exc(n);
        std::vector<double> ela(n);
        std::vector<double> tot(n);

        for(size_t i = 0; i<n; i++){
            for(const table& tab: Xsec_eff){
                const double sigma = tab.section[i];
                eff[i] += sigma;
                tot[i] += sigma;
            }
            for(const table& tab: Xsec_att){
                const double sigma = tab.section[i];
                att[i] += sigma;
                tot[i] += sigma;
            }
            for(const table& tab: Xsec_exc){
                const double sigma = tab.section[i];
                exc[i] += sigma;
                tot[i] += sigma;
            }
            for(const table& tab: Xsec_ela){
                const double sigma = tab.section[i];
                ela[i] += sigma;
                tot[i] += sigma;
            }
        }

        // Import data in files:
        std::stringstream ss;
        ss << "plots/Xsec/" << species[spec_idx] << ".csv";
        std::string filename = ss.str();

        std::ofstream file(filename);
        for(size_t i = 0; i < n; i++) {
            file << eff[i] << "," << att[i] << "," << exc[i] << "," << ela[i] << "," << tot[i] << "\n";
        }
        file.close();
    }
}