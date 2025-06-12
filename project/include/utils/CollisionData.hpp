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
    // Constructors:
    CollisionData() = default;
    CollisionData(const CrossSectionsData& Xsec, const std::vector<double>& mgas);

    // Public Method:
    void ComputeIndeces(const int n_electrons, const CrossSectionsData& Xsec,
        const std::vector<double>& E_in_eV, const std::vector<double>& mix,
        const double density, const std::vector<double>& R);

    // Getters:
    const std::vector<double>& getMass() const { return Mass; }
    const std::vector<double>& getLoss() const { return Loss; }
    const std::vector<size_t>& get_ind(std::string type) const;
    const int getCollisions() const;

private:
    
    // Class Members:

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

    // Private Methods:
    const size_t CollisionMatrix(const double R, const CrossSectionsData& Xsec,
        const std::vector<double>& mix, const double E_in_eV, const double density);
    void fill_Mass(const std::vector<size_t>& ind);
    void fill_Loss(const std::vector<size_t> & ind);
    void find_collision_indeces( const std::vector<size_t>& ind);

};

#endif // COLLISION_MATRIX_HPP