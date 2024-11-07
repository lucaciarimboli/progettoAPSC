#ifndef MONTECARLO_H
#define MONTECARLO_H

#include <string>
#include <vector>
#include <unordered_set>

class cross_sect;

enum ParticleType {ELECTRONS = 0, CATIONS, ANIONS, PARTICLES_TYPES};
struct MeanData     // Mean values for electrons
{
        double energy;
        std::array<double, 3> position;
        std::array<double, 3> sigma;
        std::array<double, 3> velocity;
        std::array<int, PARTICLES_TYPES> particles; // # of particles per each type
};


class MonteCarlo
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
    cross_sect Xsec;
    // cell array of subformula of gas species
    std::vector<std::string> gas;
    // cell array of mass of gas species (in kg)
    std::vector<double> mgas;
    // fractions of individual species in the gas as a vector
    std::vector<double> mix;
            
    // number of initial electrons used in MC calculation
    const int N0;
    // number of initial electrons used for space charge calculation
    const int n0;
    // pressure in Pascal
    double p;
    // voltage in V
    double U;
    // distance in m
    double d;
    // temperature in Kelvin
    double Temp;
            
    // gas density in m^-3
    double N;

    // E/N for homogeneous field in Td
    EN = [];
    // minimal E/N for inhomogeneous field in Td
    EN_min = [];
    // maximal E/N for inhomogeneous field in Td
    EN_max = [];
            
            
    // length in x direction
    const double Lx;
    // length in y direction
    const double Ly;
    // length in z direction
    const double Lz;
    // number of cells in x direction
    const int nx = 80;
    // number of cells in y direction
    const int ny = 90;
    // number of cells in z direction
    const int nz = 100;
    // x_vector
    std::array<double,nx+1> x;
    // y_vector
    std::array<double,ny+1> y;
    // z_vector
    std::array<double,nz+1> z;
    // x_meshgrid
    X;
    // y_meshgrid
    Y;
    // z_meshgrid
    Z;
    // [nx x ny x nz]-matrix with zeros inside and ones outside the boundary
    // boundary;
            
    // array of initial mean position of initial gaussian distributed electrons in x,y and z direction
    const std::array<double,3> pos_xyz  = {0, 0, 0};
    // array of initial broadening of initial gaussian distributed electrons in x,y and z direction
    const std::array<double,3> sigma_xyz = {0, 0, 0};
                    
    // tolerance in error of drift velcocity (default: 1%)
    static double w_err = 0.001;
    // tolerance in error of diffusion constant (default: 1%)
    static double DN_err = 0.001;
    // maximum allowed number of electorns
    static double Ne_max = 1e6;
    // number of collisions until equilibrium (default: 20e6)
    static double col_equ = 10e6;
    // number of collisions at which simulation ends (default: 20e6) 
    static double col_max = 20e6;
    // conserve (1) electron number after ionizatzion/attachment or not (0)
    bool conserve = 1;    
    // (1) isotropic, (0) non-isotropic scattering according to Vahedi et al.
    bool iso = 1;
    // energy sharing in ionizing collision
    W;
    // maximum electron energy
    E_max;
    // maximal collision frequency:
    nu_max;
    // collision counter
    counter;
    // checks end of simulation: End =1 stops the simulation
    End = 0;
    // euqilibrium time
    T_sst;
    // line number in output file:
    line = 1;
    // computation time:
    elapsedTime;
    // plot (1) or do not plot data (0)
    interactive;
            
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
    // current time-integrated velocity
    v_int;
    // current time-integrated velocity-squared
    v2_int;
    // collision indices for elastic collision
    ind_ela;
    // collision indices for excitation collision
    ind_exc;
    // collision indices for ionization collision
    ind_ion;
    // collision indices for attachment collision
    ind_att;
    // total number of all real collisions that happend
    collisions = 0;
    // column numbers of elastic collision
    col_ela;
    // column numbers of excitation collision
    col_exc;
    // column numbers of ionization collision
    col_ion;
    // column numbers of attachment collision
    col_att;
    // Mass is a vector of length(v), the entries are the masses of the gas-species undergoing elastic
    // collisions, all other entries are zero
    Mass;
    // Loss is a vector of length(v), the entries are the energy losses due to excitation collisions,
    //all other entries are zero
    Loss;
            
    // temporal mean data of electron swarm
    MeanData mean;
    // bulk transport data
    bulk;
    // flux transport data
    flux;
    // reaction rates
    rates;
    // energy data
    E;
            
    // charge density
    std::array<std::vector<double>,PARTICLES_TYPES> rho;        // DA CORREGGERE !!
    // electric potential
    phi;
    // electric field function in x-direction
    E_x;
    // electric field function in y-direction
    E_y;
    // electric field function in z-direction
    E_z;
            
    //Heat energy partition
    EnergyLossElastic    = 0;
    EnergyLossInelastic  = 0;
    EnergyLossIonization = 0;
            
    //Status
    bool converge = 0;

    // DEFAULT CONSTRUCTOR, DESTRUCTOR, COPY-CONSTRUCTOR
    // BUILD CONSTRUCTOR TO IMPORT VALUES AS IN 'MonteCarlo_singlerun'

    // Checks if sum of gas fractions is equal to 1
    // if not the case: last entry of mix will be corrected
    void checkFractionSum();
    
    // Calculates the molar weight of a substance
    std::vector<double> MolMass(std::string & substance);
    // Auxilarry methods for the "MolMass" method:
    void check_syntax(const std::string & substance, std::unordered_set<char> & characters);
    void fix_spaces(std::string & substance);
    struct Element {    // type needed in "MolMass" method
        std::string substance_char;     // element chemical symbol (e.g. "O","Ar")
        double substance_mass;          // single substance mass
        unsigned int factor = 1;        // single substance number factor (e.g. "2" in "O2")
    };
    
    // Converts mass (mgas) for the gas species from a.u. into kg
    void mass_in_kg();

    // Calculates the gas number density (in m^-3) by the ideal
    // gas law from pressure p and temperature Temp
    void gasNumberDensity();

    // Calculates absolute value of velocity and energy in eV
    std::pair<double,double> velocity2energy_in_ev(const std::vector<double> & v);

    // Calculates maximal collision rate for gas mixture (in s^-1) 
    // void maximalCollFreq();                                // USA CLASSE "cross_sect"

    // Sets initial position and velocity of electrons
    void initialParticles();

    // Makes a 3d-meshgrid and defines boundaries       // CLASS MESH ? 
    // void geometry();                                 // MIGHT BE UNUSED ??

    // Removes electrons outside the boundary 
    void surfaceInteraction();

    // Calculates electron number density (in 1/m^3)     
    void particle2density();




    
};

#endif 