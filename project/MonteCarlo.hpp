#ifndef MONTECARLO_H
#define MONTECARLO_H

#include <string>
#include <vector>
#include <unordered_set>

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

    // cross secion data created by class importLXcat
    Xsec;
    // cell array of sumformula of gas species
    std::vector<std::string> gas;
    // cell array of mass of gas species (in kg)
    std::vector<double> mgas;
    // fractions of individual species in the gas as a vector
    std::vector<double> mix;
            
    // number of initial electrons used in MC calculation
    int N0;
    // number of initial electrons used for space charge calculation
    int n0;
    // pressure in Pascal
    double p;
    // voltage in V
    U;
    // distance in m
    d;
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
    Lx;
    // length in y direction
    Ly;
    // length in z direction
    Lz;
    // number of cells in x direction
    nx = 80;
    // number of cells in y direction
    ny = 90;
    // number of cells in z direction
    nz = 100;
    // x_vector
    x;
    // y_vector
    y;
    // z_vector
    z;
    // x_meshgrid
    X;
    // y_meshgrid
    Y;
    // z_meshgrid
    Z;
    // [nx x ny x nz]-matrix with zeros inside and ones outside the boundary
    boundary;
            
    // vector of initial mean position of initial gaussian distributed electrons in x,y and z direction
    pos_xyz  = [0 0 0];
    // vector of initial broadening of initial gaussian distributed electrons in x,y and z direction
    sigma_xyz = [0 0 0];
            
            
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
    conserve = 1;    
    // (1) isotropic, (0) non-isotropic scattering according to Vahedi et al.
    iso = 1;
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
    t = [];
    // current time step dt
    dt;
    // sum of all times for all electrons:
    t_total = 0;
    // current position of electrons, cations and anions (order important)
    r = {[] [] []};
    // current velocity of electrons
    v;
    // current acceleration of electrons
    a;
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
    mean;
    // bulk transport data
    bulk;
    // flux transport data
    flux;
    // reaction rates
    rates;
    // energy data
    E;
            
    // charge density
    rho;
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
    converge = 0;

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

    // calculates the gas number density (in m^-3) by the ideal
    // gas law from pressure p and temperature Temp
    void gasNumberDensity();

    // calculates absolute value of velocity and energy in eV
    std::pair<double,double> velocity2energy_in_ev(const std::vector<double> & v);

    
};

#endif 