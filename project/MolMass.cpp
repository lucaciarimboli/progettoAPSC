#include "MolMass.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <cmath>

void MolMass::check_syntax(){
    for( auto it = substance.cbegin(); it != substance.cend(); it++){
        if( characters.find(*it) == characters.cend() ){
            std::cerr << "Wrong symbol '" << *it << "' in the formula!" << std::endl;
        }
    }
}

void MolMass::fix_spaces(){
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

void MolMass::Compute_M(){
    // "substance" is a string of the chemical formula of s substance.
    // example "Fe2(SO4)3";
    // In this case "M" would have length 1: it contains the molar mass of the substance

    // "substance" can also be a vector of substances opened by '[' and divided by space, comma or semicolon.
    // examples:
    // "[Fe2(SO4)3 CuSO4 NaOH]";
    // "[H2SO4;H2O;P;Cl2]";
    // "[C3H5(OH)3,C3H7OH,C12H22O11,NaCl]";
    // In this case "M" would have as components the molar masses of each susbtance. 

    check_syntax(); // Check "substance" syntax
    fix_spaces(); // Fix spaces within "substance"

    auto it = substance.cbegin();

    M.clear(); 
    
    if( *it == '['){    // GESTISCI NEL CONSTRUCTOR!!!

        it++;    // skip the first char '['
        std::string single_substance;

        while( it != this->substance.cend() ){
            if( *it == ' ' || *it == ',' || *it == ';' || *it == ']' ){
                MolMass single_molmass(single_substance);
                single_molmass.Compute_M();
                this->M.emplace_back(single_molmass.get_front_M());    // MolMass(single_substance) returns a vector<double> of size 1, not a double --> .front() needed
                single_substance.clear();
                while( it != this->substance.cend() && 
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

        // substance string might be either written as a name or as a chemical formula:
        if( *it != '"'){ // " indicates substance written as a name

            // A std::vector stores separately substance single elements as MolMass objects.
            // e.g. if substance="H20", elements=[(H2 data);(O data)]
            std::vector<MolMass> elements; 

            bool number_flag = 0;   // necessary for computing element.factor

            while( it != substance.cend()){
                if((*it >= 'A') && (*it <= 'Z')){
                    MolMass elem;
                    elem.append_substance(*it);
                    // elem.factor = 1;
                    elements.emplace_back(elem);
                    number_flag = 0;
                }
                else if((*it >= 'a') && (*it <= 'z')){
                    MolMass & elem = elements.back(); // elem is the last Element created
                    elem.append_substance(*it);
                    // elements.back().factor = 1;
                    number_flag = 0;
                }
                else if((*it >= '0') && (*it <= '9')){
                    MolMass & elem = elements.back();
                    if(number_flag){    // The previous char of 'substance' was a number
                        elem.factor = 10*elem.factor + *it - '0'; 
                    } else{
                        elem.factor = *it -'0';
                    }
                    number_flag = 1;
                }
                else if((*it == '+') || (*it == '-')){
                    MolMass & elem = elements.back();
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
                    // One Element accounts for the whole substance within () is added to "elements":
                                       
                    MolMass elem(bracket_substance);
                    elem.Compute_M(); 
                    elem.set_substance("$");    // Mark that "elem" is actually a substance that was inside brackets.
                    elements.emplace_back(elem); // Add brackets substance to elements
                    number_flag = 0;
                }
                it++;
            }
            
            // Update molar masses of elements from "elements_mm" (unordered_map with mm data)
            for(MolMass & elem : elements){
                if(elem.get_substance() != "$"){     // for substances in brackets, "substance_mass" is already filled
                        auto mm_iterator = this->elements_mm.find(elem.get_substance());
                        elem.set_M(mm_iterator->second);
                }
            }

            double computed_MM = 0;

            // Compute total molar mass
            for(const MolMass & elem : elements){
                computed_MM += elem.factor * elem.get_front_M();
            }
            this->M.emplace_back(computed_MM);   // MM is a vector<double> to account for the case of multi-substances input
        }

        else if(this->substance.find_first_of("(") != std::string::npos){     // substance is not a chemical formula but a name
            for( size_t i = 1; i <= this->substance.length(); i++){
                if(this->substance[i] == '('){
                    size_t j = i + 1;
                    std::string MM_str;
                    while( this->substance[j] != ')'){
                        MM_str += this->substance[j];
                        j++;
                    }
                    this->M.emplace_back(stod(MM_str));
                }
            }
        }

        else this->M.emplace_back(1);
    }
}