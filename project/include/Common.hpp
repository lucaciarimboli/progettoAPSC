#ifndef COMMON_HPP
#define COMMON_HPP

#include <vector>
#include <array>
#include <string>



namespace mc {

    //---------------------------//
    //         Constants         //
    //---------------------------//

    // electron mass
    inline constexpr double me = 9.10938291e-31;
    // electron charge
    inline constexpr double q0 = 1.60217657e-19;
    // Boltzmann constant
    inline constexpr double kB = 1.3806488e-23;
    // Avogadro constant
    inline constexpr double Na = 6.02214129e23;
    // electric constant
    inline constexpr double epsilon0 = 8.854188e-12;
    // factor 1/2 m/q to convert velocity into energy (eV)
    inline constexpr double factor = 0.5 * me / q0;

    //---------------------------//
    //          Types            //
    //---------------------------//

    /// Alias for 3D positions and velocities of particles populations.
    /// Every element is a 3 doubles array referring to x,y and z directions.
    using MATRIX =  std::vector<std::array<double,3>>;

    /**
     * @enum ParticleType
     * @brief Types of particles.
     * 
     * The last particle type does not represents a particle,
     * but it corresponds to the number of particle types.
     * It can be used in loop or when defining container sizes
     * to clarify that the indexes 0,1,2 refer to particle types.
     * (e.g. std::array<int,PARTICLES_TYPES> instead of std::array<int,3>)
     */
    enum ParticleType {
        ELECTRONS = 0,  ///< electrons.
        CATIONS,        ///< cations.
        ANIONS,         ///< anions.
        PARTICLES_TYPES ///< number of particle types.
    };

    /**
     * @enum InteractionType
     * @brief Types of interactions.
     * 
     * The last interaction type does not represents an interaction,
     * but it corresponds to the number of interaction types.
     * It can be used in loop or when defining container sizes
     * to clarify that the indexes 0,1,2 refer to particle types.
     */
    enum InteractionType {
        EFFECTIVE = 0,  ///< effective.
        IONIZATION,     ///< ionization.
        ATTACHMENT,     ///< attachment.
        EXCITATION,     ///< excitation.
        ELASTIC,        ///< elastic.
        INTERACTIONS    ///< number of interaction types.
    };
}

# endif // COMMON_HPP