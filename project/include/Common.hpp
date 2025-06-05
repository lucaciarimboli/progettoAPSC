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

    //---------------------------//
    //          Types            //
    //---------------------------//

    // MATRIX is used for 3D positions or velocities of particles populations:
    using MATRIX =  std::vector<std::array<double,3>>;

    // To improve clarity when referring to the different types of particles
    enum ParticleType {
        ELECTRONS = 0,
        CATIONS,
        ANIONS,
        PARTICLES_TYPES
    };

    // To improve clarity when referring to the different types of interactions
    enum InteractionType {
        EFFECTIVE = 0,
        IONIZATION,
        ATTACHMENT,
        EXCITATION,
        ELASTIC,
        INTERACTIONS
    };
}

# endif // COMMON_HPP