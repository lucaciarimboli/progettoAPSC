//#include "MonteCarlo.hpp"
#include "test_class.hpp"

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

void MonteCarlo::check_syntax(const std::string & substance, std::unordered_set<char> & characters){
    for( auto it = substance.cbegin(); it != substance.cend(); it++){
        if( characters.find(*it) == characters.cend() ){
            std::cerr << "Wrong symbol '" << *it << "' in the formula!" << std::endl;
        }
    }
}

void MonteCarlo::fix_spaces(std::string & substance){
    // Keep at most one space in between chars
    size_t pos = substance.find("  ");
    while (pos != std::string::npos) {
        substance.replace(pos, 2, " ");
        pos = substance.find("  ");
    }

    // Remove spaces around brackets
    pos = substance.find(" ]");
    if (pos != std::string::npos) {
        substance.replace(pos, 2, "]");
    }
    pos = substance.find("[ ");
    if (pos != std::string::npos) {
        substance.replace(pos, 2, "[");
    }
}

std::vector<double> MonteCarlo::MolMass(std::string & substance){
    // "substance" is a string of the chemical formula of s substance.
    // example:	MM = MolMass('Fe2(SO4)3');

    // "substance" can also be a vector of substances opened by '[' and divided by space, comma or semicolon.
    // examples:
    // MM = MolMass('[Fe2(SO4)3 CuSO4 NaOH]');
    // MM = MolMass('[H2SO4;H2O;P;Cl2]');
    // MM = MolMass('[C3H5(OH)3,C3H7OH,C12H22O11,NaCl]');

    // To distinguish charched substances the symbols '+' and '-' can be used.
    // exampels:
    // MM = MolMass('Fe2+')  --->  MM = 55.8470		      (it means one mol of Fe2+)
    // MM = MolMass('Fe3+')  --->  MM = 55.8470		      (it means one mol of Fe3+)
    // but	MM = MolMass('Fe2')   --->  MM = 111.6940	  (it means two moles of Fe)

    std::unordered_set<char> characters{'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 
                                        'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 
                                        'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 
                                        'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', 
                                        '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '(', ')', '[', 
                                        ']', ';', ',', ' ', '+', '-', '.'};

    check_syntax(substance,characters); // Check "substance" syntax
    fix_spaces(substance); // Fix spaces within "substance"

    std::vector<double> MM; // Vector with molar masses to be returned;

    auto it = substance.cbegin();

    if( *it == '['){    // Vector of substances

        it++;    // skip the first char '['
        std::string single_substance;

        while( it != substance.cend() ){
            if( *it == ' ' || *it == ',' || *it == ';' || *it == ']' ){
                MM.emplace_back(MolMass(single_substance).front());    // MolMass(single_substance) returns a vector<double> of size 1, not a double --> .front() needed
                single_substance.clear();
                while( it != substance.cend() && 
                    (*it == ' ' || *it == ',' || *it == ';' || *it == ']') ){
                    it++;
                }
            } else if( characters.find(*it) != characters.cend() ){
                single_substance += *it;
                it++;
            }
        }
    }

    else{       // Single substance
        std::unordered_map<std::string, double> elements_mm = {
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

        // substance string might be either written as a name or as a chemical formula:
        if( *it != '"'){ // " indicates substance written as a name

            // A std::vector stores separately substance elements data.
            // e.g. if substance="H20", elements=[(H2 data);(O data)]
            std::vector<Element> elements_vect; 

            bool number_flag = 0;   // necessary for computing element.factor

            while( it != substance.cend()){
                if((*it >= 'A') && (*it <= 'Z')){
                    Element elem;
                    elem.substance_char += *it;
                    // elem.factor = 1;
                    elements_vect.emplace_back(elem);
                    number_flag = 0;
                }
                else if((*it >= 'a') && (*it <= 'z')){
                    Element & elem = elements_vect.back(); // elem is the last Element created
                    elem.substance_char += *it;
                    // (elements_vect.back()).factor = 1;
                    number_flag = 0;
                }
                else if((*it >= '0') && (*it <= '9')){
                    Element & elem = elements_vect.back();
                    if(number_flag){    // The previous char of 'substance' was a number
                        elem.factor = 10*elem.factor + *it - '0'; 
                    } else{
                        elem.factor = *it -'0';
                    }
                    number_flag = 1;
                }
                else if((*it == '+') || (*it == '-')){
                    Element & elem = elements_vect.back();
                    elem.factor = 1;
                }
                else if(*it == '('){
                    unsigned int bracket_level = 1;     // keeps track of the level of nested parentheses
                    std::string bracket_substance;  // string to store the substance between "()""

                    while( bracket_level > 0){
                        it++;
                        if(*it == '(') bracket_level++;
                        else if (*it == ')') bracket_level--;
                        
                        if( bracket_level > 0) bracket_substance += *it;
                    }
                    // One Element accounts for the whole substance within () is added to "elements_vector":
                                       
                    Element elem; 
                    elem.substance_mass = MolMass(bracket_substance).front(); // Recursive call to compute MM of substance inside brackets
                    elem.substance_char = "(";    // Mark that "elem" is actually a substance that was inside brackets.
                    elements_vect.emplace_back(elem); // Add brackets substance to elements_vect
                    number_flag = 0;
                }
                it++;
            }
            
            // Update molar masses of elements from "elements_mm" (unordered_map with mm data)
            for(Element & elem : elements_vect){
                if(elem.substance_char != "("){     // for substances in brackets, "substance_mass" is already filled
                        auto mm_iterator = elements_mm.find(elem.substance_char);
                        elem.substance_mass = mm_iterator->second;
                }
            }

            double computed_MM = 0;

            // Compute total molar mass
            for(const Element & elem : elements_vect){
                computed_MM += elem.factor * elem.substance_mass;
            }
            MM.emplace_back(computed_MM);   // MM is a vector<double> to account for the case of multi-substances input
        }

        else if(substance.find_first_of("(") != std::string::npos){     // substance is not a chemical formula but a name
            for( size_t i = 1; i <= substance.length(); i++){
                if(substance[i] == '('){
                    size_t j = i + 1;
                    std::string MM_str;
                    while( substance[j] != ')'){
                        MM_str += substance[j];
                        j++;
                    }
                    MM.emplace_back(stod(MM_str));
                }
            }
        }

        else MM.emplace_back(1);
    }
    return MM;
}

// MolMass returns a std::vector<double> since e.g. with an input
// "[H2O,CO2]" as one string it returns a vector of two doubles with the
// molar mases of, respectively, "H2O" and "CO2".
// For this reason "mass_in_kg" has two for loops, even if most likely in the 
// application only one substance per string would be in "gas". 

void MonteCarlo::mass_in_kg(){
    // converts mass (mgas) for the gas species from a.u. into kg.
    this->mgas.clear(); // In case the method is called more than once by mistake
    for(const std::string & g : this->gas){
        for( const double & mm : this->MolMass(g / ( this->Na * 1000))){
            this->mgas.emplace_back(mm);
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