#ifndef MEANDATA_H
#define MEANDATA_H

#include <array>
#include <vector>
#include <iostream>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "Common.hpp"

/**
 * @class MeanData
 * @brief Computes and stores instantaous mean values of an electron swarm.
 *
 * This class stores the mean characteristics (position, velocity, energy, variance)
 * of a particle population, including the number of particles per type.
 * Types of particles are: electrons, cations, anions.
 * The mean values refers to a specific instant in time, the averaging is performed
 * among the positions, velocities and energies of the electron swarm at a given time.
 */
class MeanData
{
    public:

    /**
     * @brief Constructor for initial condition.
     * @param r Mean position vector.
     * @param s Variance of position.
     * @param ne Number of electrons.
     * 
     * Constructs MeanData from known mean position, variance and number of electrons.
     */
    MeanData(const std::array<double,3>& r, const std::array<double,3>& s, const int& ne);
    /**
     * @brief Constructor for t>0.
     * @param p Number of particles per each type
     * @param r Matrix of positions of every electron.
     * @param v Matrix of velocities of every electron.
     * @param E_in_eV Vector of kinetic energies of every electron.
     * 
     * Constructs MeanData by averaging the electron swarm data provided.
     */
    MeanData(const std::array<int,mc::PARTICLES_TYPES>& p, const mc::MATRIX& r, const mc::MATRIX& v, const std::vector<double>& E_in_eV);

    /**
     * @brief Updates the number of particles by adding new counts.
     * @param p Number of new particles per type to add.
     */
    void add_new_particles(const std::array<int,mc::PARTICLES_TYPES>& p);

    /**
     * @brief Sets the number of particles per type.
     * @param p Number of particles per type.
    */
    void set_particles(const std::array<int,mc::PARTICLES_TYPES>& p){ particles = p; };

    /**
     * @brief Sets the electrons center of mass.
     * @param pos Center of mass coordinates
    */
    void set_position(const std::array<double,3>& pos){ position = pos; };

    /**
     * @brief Sets the electrons variance in position.
     * @param v Variance of electron positions.
    */
    void set_variance(const std::array<double,3>& v){ var = v; };


    // Getters:
    /**
     * @brief Gets the number of particles per type.
     * @return Number of particles array.
     */
    const std::array<int,mc::PARTICLES_TYPES> & get_particles() const{ return particles; }
    /**
     * @brief Gets the center of mass of the electrons.
     * @return Electrons center of mass array.
     */
    const  std::array<double,3>& get_position() const{ return position; }
    /**
     * @brief Gets the variance of the electrons positions.
     * @return Electrons position variance array.
     */
    const std::array<double,3>& get_variance() const{ return var; }
    /**
     * @brief Gets the mean velocity of the electrons.
     * @return Electrons mean velocity array.
     */
    const std::array<double,3>& get_velocity() const{ return velocity; }
    /**
     * @brief Gets the mean kinetic energy of the electrons.
     * @return Electrons mean kinetic energy.
     */
    double get_energy() const{ return energy; }

    private:

    // Class Members:
    std::array<int,mc::PARTICLES_TYPES> particles; ///< Number of particles per each type
    std::array<double, 3> position;                ///< mean position
    std::array<double, 3> var;                     ///< position variance
    std::array<double, 3> velocity;                ///< mean velocity
    double energy;                                 ///< mean energy [eV]
};

#endif