#ifndef COLLISION_MATRIX_H
#define COLLISION_MATRIX_H

#include <vector>
#include <algorithm>
#include <unordered_set>
#include <random>

#include "utils/CrossSectionsData.hpp"
#include "Common.hpp"

class CollisionData {

public:
    // Constructor
    CollisionData() = default;
    CollisionData(const CrossSectionsData& Xsec, const std::vector<double>& mgas);

    // Build collision matrix and compute indeces:
    void ComputeIndeces(const int n_electrons, const CrossSectionsData& Xsec,
        const std::vector<double>& E_in_eV, const std::vector<double>& mix,
        const double density, const std::vector<double>& R);

    // Getters:
    //const double getC(const int i, const int j, const int n_articles) const { return C[i*n_particles + j]; }
    const std::vector<double>& getMass() const { return Mass; }
    const std::vector<double>& getLoss() const { return Loss; }
    const std::vector<size_t>& get_ind(std::string type) const {
        if (type == "ELASTIC")         return ind_ela;
        else if (type == "EXCITATION") return ind_exc;
        else if (type == "IONIZATION") return ind_ion;
        else if (type == "ATTACHMENT") return ind_att;
        else throw std::invalid_argument("Invalid collision type");
    }
    // Count and return the number of real collisions happened:
    const int getCollisions() const {
        return ind_ela.size() + ind_exc.size() + ind_ion.size() + ind_att.size();
    }

private:

    // Vectors with fixed indeces defined by the constructor:
    std::vector<double> mass;             // gas specie mass for each reaction
    std::vector<double> loss;             // energy loss for each reaction

    std::unordered_set<size_t> col_ela;           // react indeces for elastic collision
    std::unordered_set<size_t> col_exc;           // react indeces for excitation collision
    std::unordered_set<size_t> col_ion;            // react indeces for ionization collision
    std::unordered_set<size_t> col_att;            // react indeces for attachment collision

    // Vectors to be updated at each time step:
    std::vector<double> Mass;              // Mass vector
    std::vector<double> Loss;              // Loss vector    
    
    std::vector<size_t> ind_ela;           // Collision indices for elastic collision
    std::vector<size_t> ind_exc;           // Collision indices for excitation collision
    std::vector<size_t> ind_ion;           // Collision indices for ionization collision
    std::vector<size_t> ind_att;           // Collision indices for attachment collision

    // Build collision matrix and update current index based on the random number:
    const size_t CollisionMatrix(const double R, const CrossSectionsData& Xsec,
        const std::vector<double>& mix, const double E_in_eV, const double density);

    // Fill the Mass vector with the mass of the gas species:
    void fill_Mass(const std::vector<size_t>& ind);
    // Fill the Loss vector with the energy loss of the gas species:
    void fill_Loss(const std::vector<size_t> & ind);
    // Find the collision indeces for each type of collision:
    void find_collision_indeces( const std::vector<size_t>& ind);

};

#endif // COLLISION_MATRIX_HPP