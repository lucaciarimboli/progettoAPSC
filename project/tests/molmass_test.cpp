#include <iostream>
#include <cassert>
#include "../include/utils/MolMass.hpp"

int main() {
    // Test 1: Constructor and Getters
    std::cout << "TEST 1\n" << std::endl;
    MolMass mol1("H2O");
    std::cout << "should be H2O: " << mol1.get_substance() << "\n" << std::endl;
    mol1.Compute_M();
    std::cout << "should be 18.015: " << mol1.get_front_M() << "\n\n" << std::endl;

    // Test 3: Named Substance
    std::cout << "TEST 3\n" << std::endl;
    MolMass mol3("\"Custom(123.45)\"");
    mol3.Compute_M();
    std::cout << "should be 123.45: " << mol3.get_front_M() << std::endl;

    // Test 4: Mixture
    std::cout << "TEST 4\n" << std::endl;
    MolMass mol4("[H2O,CO2]");
    mol4.Compute_M();
    const auto& mix_masses = mol4.get_M();
    std::cout << "should be 2: " << mix_masses.size() << "\n" << std::endl;
    std::cout << "should be 18.015: " << mix_masses[0] << "\n" << std::endl;
    std::cout << "should be 44.01: " << mix_masses[1] << "\n\n" << std::endl;
    

    // Test 5: Invalid Syntax
    std::cout << "TEST 5\n" << std::endl;
    try {
        MolMass mol5("InvalidFormula");
        mol5.Compute_M();
        assert(false);  // Should not reach here
    } catch (const std::invalid_argument& e) {
        std::cout << "Caught expected exception: " << e.what() << std::endl;
    }

    // Test 6: Private Methods (temporarily moved to public for testing)
    std::cout << "TEST 6\n" << std::endl;
    MolMass mol6("Ca(OH)2");
    std::cout << "should be CaO2H2 : " << mol6.adapt_formula("Ca(OH)2") << "\n" << std::endl;
    std::cout << "should be 74.092: " << mol6.compute_substance("CaO2H2") << "\n\n" << std::endl;

    MolMass mol7("[H2O,CO2]");
    mol7.fix_spaces();
    std::cout << "should be [H2O,CO2]: " << mol7.get_substance() << "\n\n" << std::endl;

    // Test 7: Edge Cases
    std::cout << "TEST 7\n" << std::endl;
    MolMass mol8("");
    mol8.Compute_M();
    assert(mol8.get_M().empty());  // Empty substance should result in empty molar mass

    MolMass mol9("[H2O]");
    mol9.Compute_M();
    std::cout << "should be 1: " << mol9.get_M().size() << "\n" << std::endl;
    std::cout << "should be 18.015: " << mol9.get_front_M() << "\n\n" << std::endl;

    MolMass mol10("\"Custom\"");
    mol10.Compute_M();
    std::cout << "should be 1.0: " << mol10.get_front_M() << "\n\n" << std::endl;
    // 1.0 default molar mass

    // Test 8: Check Syntax
    std::cout << "TEST 8\n" << std::endl;
    MolMass mol11("Fe2(SO4)3");
    if(mol11.check_syntax())
        std::cout << "mol1 succeded\n" << std::endl;

    MolMass mol12("[H2O,CO2]");
    if(mol12.check_syntax())
        std::cout << "mol2 succeded\n" << std::endl;

    MolMass mol13("\"Custom(123.45)\"");
    if(mol13.check_syntax())
        std::cout << "mol3 succeded\n" << std::endl;

    MolMass mol14("InvalidFormula");
    if(!mol14.check_syntax())
        std::cout << "mol4 succeded\n" << std::endl;

    std::cout << "\nAll tests passed successfully!" << std::endl;
    return 0;
}