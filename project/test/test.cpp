#include "test_class.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <cmath>

int main(){
    
    MonteCarlo mc;
    mc.substance = "[PO4 , H2O , Ar3]";
    //mc.substance = "Ca3(PO4)2"; // returns 310.177 CORRECT! 
    //mc.substance = "(PO4)2"; // segmentation fault

    auto molmass = mc.MolMass(mc.substance);

    std::cout << "Substance = " << mc.substance << "\n Molar Mass = " << std::endl;
    for(const double & mm : molmass){
        std::cout << mm << ", " << std::endl;
    }

    return 0;
}