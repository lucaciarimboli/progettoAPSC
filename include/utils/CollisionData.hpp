#ifndef COLLISION_MATRIX_H
#define COLLISION_MATRIX_H

#include <vector>
#include <algorithm>
#include <unordered_set>
#include <random>

#include "utils/CrossSectionsData.hpp"
#include "Common.hpp"

// #include <omp.h>

/**
 * @brief Class to compute and store collision data for a Monte Carlo electron transport simulation.
 * 
 * This class handles the processing of cross section data and electron transport data
 * to define a collision matrix and, provided a random number for each electron of the swarm, 
 * assign to every electron a collision event with a backround gas particle
 * The types of collisions are: elastic, excitation, ionization, attachment.
 * Also, an electron might be assigned to not collide with gas particles.
 */
class CollisionData {

public:

    /**
     * @brief Default constructor.
     */
    CollisionData() = default;

    /**
     * @brief Constructor.
     * 
     * @param Xsec Cross section data for all species and reactions.
     * @param mgas Molecular mass of each gas species in the mixture.
     */
    CollisionData(const CrossSectionsData& Xsec, const std::vector<double>& mgas);

    /**
     * @brief Compute the collision indices for a given electron population.
     * 
     * @param n_electrons Number of electrons in the simulation.
     * @param Xsec Cross section data.
     * @param E_in_eV Vector of electron kinetic energies.
     * @param v_abs Vector of absolute velocities for each electron.
     * @param mix Vector of gas mixture fractions.
     * @param density Background gas number density.
     * @param R Vector of random numbers for collision selection (size = n_electrons).
     * 
     * Given a vector containing one random number from a uniform distribution in [0,1] for each electron,
     * given the required cross sections and electron transport data, assigns to every electron
     * a collision and stores the electron indices in apposite vectors based on the type of
     * collision they have assigned to undergo.
     * 
     */
    void AssignCollisions(const int& n_electrons, const CrossSectionsData& Xsec,
        const std::vector<double>& E_in_eV, const std::vector<double>& v_abs, const std::vector<double>& mix,
        const double& density, const std::vector<double>& R);

    // Getters:
    /**
     *  @brief Get the mass vector.
     *  @return Mass vector.
     */
    const std::vector<double>& getMass() const { return Mass; }
    /**
     *  @brief Get the loss vector.
     *  @return Loss vector.
     */
    const std::vector<double>& getLoss() const { return Loss; }
    /**
     * @brief Get the vector of indices for the specified collision type.
     * 
     * @param type Collision type (elastic, excitation, etc.).
     * @return Vector of indices of electrons undergoing the specified collision.
     */
    const std::vector<size_t>& get_ind(mc::InteractionType type) const;
    
    /**
     * @brief Get the total number of collisions that occurred in the current timestep.
     * @return Total number of collisions.
     */
    const int getCollisions() const;

private:
    
    // Class Members:

    // Vectors with fixed indeces defined by the constructor:
    std::vector<double> mass;             ///< gas specie mass for each reaction
    std::vector<double> loss;             ///< energy loss for each reaction

    std::unordered_set<size_t> col_ela;           ///< react indeces for elastic collision
    std::unordered_set<size_t> col_exc;           ///< react indeces for excitation collision
    std::unordered_set<size_t> col_ion;           ///< react indeces for ionization collision
    std::unordered_set<size_t> col_att;           ///< react indeces for attachment collision

    // Vectors to be updated at each time step:
    std::vector<double> Mass;              ///< Mass vector
    std::vector<double> Loss;              ///< Loss vector    
    
    std::vector<size_t> ind_ela;           ///< Collision indices for elastic collision
    std::vector<size_t> ind_exc;           ///< Collision indices for excitation collision
    std::vector<size_t> ind_ion;           ///< Collision indices for ionization collision
    std::vector<size_t> ind_att;           ///< Collision indices for attachment collision

    // Private Methods:
    /**
     * @brief Compute the reaction index for a given electron.
     * 
     * @param R Random number to sample the collision type.
     * @param XS_data Vector of cross section tables for each reaction.
     * @param mix Gas mixture fractions.
     * @param t Pre-computed slope for cross section interpolation.
     * @param k Pre-computed index of the energy bin for cross section interpolation.
     * @param factor Pre computed factor v_abs * density / nu_max.
     * 
     * @return Index of the reaction assigned to the electron.
     */
    const size_t CollisionMatrix(const double& R, const std::vector<table>& XS_data, const std::vector<double>& mix, const double& t, const size_t& k, const double& factor);

    /**
     * @brief Fill the Mass vector based on the reactions indexed in `ind`.
     * 
     * @param ind Indices of the selected reactions.
     */
    void fill_Mass(const std::vector<size_t>& ind);

    /**
     * @brief Fill the Loss vector based on the reactions indexed in `ind`.
     * 
     * @param ind Indices of the selected reactions.
     */
    void fill_Loss(const std::vector<size_t> & ind);
    
    /**
     * @brief Classify the selected reaction indices into elastic, excitation, ionization, and attachment.
     * 
     * @param ind Indices of the selected reactions.
     */
    void find_collision_indeces( const std::vector<size_t>& ind);

};

#endif // COLLISION_MATRIX_HPP