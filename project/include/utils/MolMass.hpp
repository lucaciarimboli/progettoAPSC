#ifndef MOLMASS_H
#define MOLMASS_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <cmath>
#include <regex>

class MolMass
{
    public:
    // Constructors:
    MolMass()=default;
    MolMass(const std::string & s): substance(s) {};

    // Public Method:
    void Compute_M();

    // Setters:
    void set_M(const double & m){
        M.clear();
        M.emplace_back(m);
    }
    void set_substance(const std::string & s){ substance = s; }

    // Getters:
    const std::vector<double> & get_M() const{ return M; } 
    const double get_front_M() const{ return M.front(); }
    const std::string & get_substance() const{ return substance; }

    
    private:

    // Class Members:
    std::string substance;   // Substance
    std::vector<double> M;   // Molar masses in atomic units
    
    //Molar masses of each element
    const std::unordered_map<std::string, double> elements_mm = {
            {"H", 1.00794}, {"He", 4.002602}, {"Li", 6.941}, {"Be", 9.012182}, {"B", 10.811}, {"C", 12.011},
            {"N", 14.00674}, {"O", 15.9994}, {"F", 18.9984032}, {"Ne", 20.1797}, {"Na", 22.989768}, {"Mg", 24.305},
            {"Al", 26.981539}, {"Si", 28.0855}, {"P", 30.973762}, {"S", 32.066}, {"Cl", 35.4527}, {"Ar", 39.948},
            {"K", 39.0983}, {"Ca", 40.078}, {"Sc", 44.95591}, {"Ti", 47.88}, {"V", 50.9415}, {"Cr", 51.9961},
            {"Mn", 54.93805}, {"Fe", 55.847}, {"Co", 58.9332}, {"Ni", 58.69}, {"Cu", 63.546}, {"Zn", 65.39},
            {"Ga", 69.723}, {"Ge", 72.61}, {"As", 74.92159}, {"Se", 78.96}, {"Br", 79.904}, {"Kr", 83.8},
            {"Rb", 85.4678}, {"Sr", 87.62}, {"Y", 88.90585}, {"Zr", 91.224}, {"Nb", 92.90638}, {"Mo", 95.94},
            {"Tc", 98.9063}, {"Ru", 101.07}, {"Rh", 102.9055}, {"Pd", 106.42}, {"Ag", 107.8682}, {"Cd", 112.411},
            {"In", 114.82}, {"Sn", 118.71}, {"Sb", 121.75}, {"Te", 127.6}, {"I", 126.90447}, {"Xe", 131.29},
            {"Cs", 132.90543}, {"Ba", 137.327}, {"La", 138.9055}, {"Ce", 140.115}, {"Pr", 140.90765}, {"Nd", 144.24},
            {"Pm", 146.9151}, {"Sm", 150.36}, {"Eu", 151.965}, {"Gd", 157.25}, {"Tb", 158.92534}, {"Dy", 162.5},
            {"Ho", 164.93032}, {"Er", 167.26}, {"Tm", 168.93421}, {"Yb", 173.04}, {"Lu", 174.967}, {"Hf", 178.49},
            {"Ta", 180.9479}, {"W", 183.85}, {"Re", 186.207}, {"Os", 190.2}, {"Ir", 192.22}, {"Pt", 195.08},
            {"Au", 196.96654}, {"Hg", 200.59}, {"Tl", 204.3833}, {"Pb", 207.2}, {"Bi", 208.98037}, {"Po", 208.9824},
            {"At", 209.9871}, {"Rn", 222.0176}, {"Fr", 223.0197}, {"Ra", 226.0254}, {"Ac", 227.0278}, {"Th", 232.0381},
            {"Pa", 231.0359}, {"U", 238.0289}, {"Np", 237.0482}, {"Pu", 244.0642}, {"Am", 243.0614}, {"Cm", 247.0703},
            {"Bk", 247.0703}, {"Cf", 251.0796}, {"Es", 252.0829}, {"Fm", 257.0951}, {"Md", 258.0986}, {"No", 259.1009},
            {"Lr", 262.1138}, {"Rf", 263.1182}, {"Db", 262.1229}, {"Sg", 263.1182}, {"Bh", 262.1138}, {"Hs", 262.1229},
            {"Mt", 262.1229},
            {"Nn", 1.0}  // Nn = Not named
    };


    // Private Methods:
    bool check_syntax() const;
    void fix_spaces();
    std::string adapt_formula(const std::string& formula);
    double compute_substance(const std::string& formula);
    void compute_mixture();
    void compute_named_substance();
};

#endif 