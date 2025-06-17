//#include "utils/MolMass.hpp"
#include "../../include/utils/MolMass.hpp"

void MolMass::fix_spaces() {
    // The string "substance" might be a string of the chemical formula of s substance or a vector of substances
    // opened by '[', closed by ']', and divided by space, comma or semicolon.
    // This method uniforms the syntax of the string "substance" so that:
    //    - There are no spaces
    //    - In case of a mix, the separator between diffent species is ",", the opening and closure are '[', ']'.

    // Remove any space at the beginning or at the end of the string:
    substance = std::regex_replace(substance, std::regex(R"(^\s+|\s+$)"), "");

    // If the string represents a mix, remove every space and set "," as separator between different species:
    if (substance.front() == '[') {
        if(substance.back() != ']') throw std::invalid_argument("Invalid chemical formula in MolMass");

        std::string content = substance.substr(1,substance.length()-2);
        content = std::regex_replace(content, std::regex(R"(\s*[,;\s]+\s*)"), ",");
        content = std::regex_replace(content, std::regex(R"(^,+|,+$)"), "");
        substance = '[' + content + ']';
    }

    // Remove all remaining spaces:
    substance = std::regex_replace(substance, std::regex(R"(\s)"), "");
}

bool MolMass::check_syntax() const {

    // Case "substance" is a single substance:
    std::regex sub_reg(R"(^([A-Z][a-z]?\d*|\([A-Z][a-z]?\d*([A-Z][a-z]?\d*)*\)\d*)+$)");
    
    // Case "substance" is a mixture:
    std::regex mix_reg(R"(^\[[A-Z][\w(),;\s]+\]$)");
    
    // Case "substance" is a named formula:
    std::regex named_reg(R"(^".*"(\(\d+\.?\d*\))?$)");
    
    // If "substance" format does not match any of the allowed one, return false:
    return std::regex_match(substance, sub_reg) ||
           std::regex_match(substance, mix_reg) ||
           std::regex_match(substance, named_reg);
}

void MolMass::Compute_M() {
    // Computes "M" in atomic units from "substance" 

    // "substance" is a string of the chemical formula of s substance.
    // example "Fe2(SO4)3";
    // In this case "M" would have length 1: it contains the molar mass of the substance

    // "substance" can also be a vector of substances opened by '[' and divided by space, comma or semicolon.
    // examples:
    // "[Fe2(SO4)3 CuSO4 NaOH]";
    // "[H2SO4;H2O;P;Cl2]";
    // "[C3H5(OH)3,C3H7OH,C12H22O11,NaCl]";
    // In this case "M" would have as components the molar masses of each susbtance. 

    // Finally, "substance" could explicitely contain the molar mass value in the form: "Name(123.45)"
    
    M.clear();

    // Check/fix the syntax of the substance:
    if (substance.empty()) return;
    fix_spaces();
    if(!check_syntax()) throw std::invalid_argument("Invalid chemical formula in MolMass");
    
    // Handle mixture notation [sub1,sub2,...]
    if (substance.front() == '[') {
        compute_mixture();
        return;
    }

    // Handle named substance with explicit mass "Name(123.45)"
    else if (substance.front() == '"') {
        compute_named_substance();
        return;
    }
    
    // Handle chemical formula
    M.emplace_back(compute_substance(substance));
}

double MolMass::compute_substance(const std::string& formula) {
    double total_mass = 0.0;
    
    // Extend the formula by removing brackets to simplify the computation:
    std::string adapted_formula = adapt_formula(formula);
    
    // Parse elements with their counts
    std::regex reg(R"(([A-Z][a-z]?)(\d*))");
    std::sregex_iterator iter(adapted_formula.begin(), adapted_formula.end(), reg);
    std::sregex_iterator end_iter;
    
    while(iter != end_iter) {
        // Extract element string:
        std::string element = (*iter)[1].str();

        // Extract number if any:
        std::string mult_str = (*iter)[2].str();
        int multiplier = mult_str.empty() ? 1 : std::stoi(mult_str);
    
        auto it = elements_mm.find(element);
        if (it != elements_mm.end()) {
            total_mass += multiplier * it->second;
        }

        iter++;
    }
    
    return total_mass;
}

std::string MolMass::adapt_formula(const std::string& formula) {
    // Adapt the formula to allow easier molar mass computation.
    // This is done by "expanding" the brackets (if any).
    // e.g. "Ca(OH)2" --> "CaO2H2"

    std::string result = formula;
    std::regex reg(R"(\(([^()]+)\)(\d*))");
    std::smatch match;
    
    while (std::regex_search(result, match, reg)) {
        std::string content = match[1].str();
        int multiplier1 = match[2].str().empty() ? 1 : std::stoi(match[2].str());
        
        std::string expanded = "";
        std::regex element_in_brackets(R"(([A-Z][a-z]?)(\d*))");
        std::sregex_iterator iter(content.cbegin(), content.cend(), element_in_brackets);
        std::sregex_iterator end_iter;

        while(iter != end_iter) {
            std::string element = (*iter)[1].str();
            std::string mult_str = (*iter)[2].str();
            int multiplier2 = mult_str.empty() ? 1 : std::stoi(mult_str);
            expanded += element + std::to_string(multiplier2 * multiplier1);
            iter++;
        }
        
        result = std::regex_replace(result, reg, expanded, std::regex_constants::format_first_only);
    }
    
    return result;
}

void MolMass::compute_mixture() {
    // Fill the vector "M" with the molar masses of the corrisponding substances in "substance".
    // The string "substance" is of the kind: "[sub1,sub2,...,subn]" without spaces.
    std::regex reg(R"(\[([^\]]+)\])");
    std::smatch match;
    if (!std::regex_match(substance, match, reg)) return;
    std::string content = match[1].str();

    std::regex separator(R"(,)");
    std::sregex_token_iterator iter(content.cbegin(), content.cend(), separator, -1);
    std::sregex_token_iterator end_iter;
    
    while(iter != end_iter) {
        std::string formula = iter->str();
        if (!formula.empty()) M.emplace_back(compute_substance(formula));
        iter++;
    }
}

void MolMass::compute_named_substance() {
    std::regex reg(R"(\(([0-9.]+)\))");
    std::smatch match;
    
    if (std::regex_search(substance, match, reg)) {
        M.emplace_back(std::stod(match[1].str()));
    } else {
        M.emplace_back(1.0);
    }
}