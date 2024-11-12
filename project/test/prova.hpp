#ifndef PROVA_H
#define PROVA_H

#include <string>
#include <vector>
#include <unordered_set>

// class cross_sect;

enum ParticleType {ELECTRONS = 0, CATIONS, ANIONS, PARTICLES_TYPES};
struct MeanData     // Mean values for electrons
{
        double energy;
        std::array<double, 3> position;
        std::array<double, 3> sigma;
        std::array<double, 3> velocity;
        std::array<int, PARTICLES_TYPES> particles; // # of particles per each type
};


class PROVA
{
    public:

    // electron mass
    static constexpr double me = 9.10938291e-31;
    // electron charge
    static constexpr double q0 = 1.60217657e-19;
    // Boltzmann constant
    static constexpr double kB = 1.3806488e-23;
    // Avogadro constant
    static constexpr double Na = 6.02214129e23;
    // electric constant
    static constexpr double epsilon0 = 8.854188e-12;

    // cross secion data created by class cross_sect
    // cross_sect Xsec;
    // cell array of subformula of gas species
    std::vector<std::string> gas;
    // cell array of mass of gas species (in kg)
    std::vector<double> mgas;
    // fractions of individual species in the gas as a vector
    std::vector<double> mix;
            
    // number of initial electrons used in MC calculation
    // const int N0;
    // number of initial electrons used for space charge calculation
    // const int n0;
    // pressure in Pascal
    double p;
    // voltage in V
    // double U;
    // distance in m
    // double d;
    // temperature in Kelvin
    double Temp;
            
    // gas density in m^-3
    double N;

    // E/N for homogeneous field in Td
    //EN = [];
    // minimal E/N for inhomogeneous field in Td
    //EN_min = [];
    // maximal E/N for inhomogeneous field in Td
    //EN_max = [];
            
    // array of initial mean position of initial gaussian distributed electrons in x,y and z direction
    const std::array<double,3> pos_xyz  = {0, 0, 0};
    // array of initial broadening of initial gaussian distributed electrons in x,y and z direction
    const std::array<double,3> sigma_xyz = {0, 0, 0};
                    
    // tolerance in error of drift velcocity (default: 1%)
    static constexpr double w_err = 0.001;
    // tolerance in error of diffusion constant (default: 1%)
    static constexpr double DN_err = 0.001;
    // maximum allowed number of electorns
    static constexpr double Ne_max = 1e6;
    // number of collisions until equilibrium (default: 20e6)
    static constexpr double col_equ = 10e6;
    // number of collisions at which simulation ends (default: 20e6) 
    static constexpr double col_max = 20e6;
    // conserve (1) electron number after ionizatzion/attachment or not (0)
    bool conserve = 1;    
    // (1) isotropic, (0) non-isotropic scattering according to Vahedi et al.
    bool iso = 1;

            
    // current time
    std::vector<double> t;
    // current time step dt
    double dt;
    // sum of all times for all electrons:
    double t_total = 0;
    // current position of electrons, cations and anions (order important)
    // r has 3 vectors: one for each type of particles.
    // e.g. vector for e contains the position (array<double,3>) of every electron;
    // x-coordinate of i-esim electron can be accessed through: r[ELECTRONS][i][0].
    typedef std::vector<std::array<double,3>> COMPONENTS;
    std::array<COMPONENTS,PARTICLES_TYPES> r;
    // current velocity of electrons
    COMPONENTS v;
    // current acceleration of electrons
    COMPONENTS a;
 
    // charge density
    std::array<std::vector<double>,PARTICLES_TYPES> rho;        // DA CORREGGERE !!
    // electric potential
    //phi;
    // electric field function in x-direction
    //E_x;
    // electric field function in y-direction
    //E_y;
    // electric field function in z-direction
    //E_z;

    //Status
    bool converge = 0;

    // Converts mass (mgas) for the gas species from a.u. into kg
    void mass_in_kg();

    // Calculates the gas number density (in m^-3) by the ideal
    // gas law from pressure p and temperature Temp
    void gasNumberDensity();

    // Calculates absolute value of velocity and energy in eV
    std::pair<double,double> velocity2energy_in_ev(const std::vector<double> & v);  
};

#endif 