#ifndef MONTECARLO_H
#define MONTECARLO_H

#include <string>
#include <vector>
#include <array>
#include <unordered_set>

class cross_sect;
class MeanData;
class EnergyData;
class FluxData;
class BulkData;
class RateDataCount;
class RateDataConv;

enum ParticleType {ELECTRONS = 0, CATIONS, ANIONS, PARTICLES_TYPES};
// enum Coordinates{ X = 0, Y, Z, COORDINATES}
typedef std::vector<std::array<double,3>> MATRIX;

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
    std::vector<cross_sect> Xsec;
    // cell array of subformula of gas species
    std::vector<std::string> gas;
    // cell array of mass of gas species (in kg)
    std::vector<double> mgas;
    // fractions of individual species in the gas as a vector
    std::vector<double> mix;
            
    // number of initial electrons used in MC calculation
    const unsigned  N0;
    // number of initial electrons used for space charge calculation
    const unsigned  n0;
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

    // E/N for homogeneous field [Td] --> constant, uniform, user defined.
    double EN;
    // minimal E/N for inhomogeneous field in Td
    // EN_min = [];
    // maximal E/N for inhomogeneous field in Td
    // EN_max = [];
            
    // MESH CLASS --> not needed if E/N is constant+uniform.     
    // length in x direction
    // const double Lx;
    // length in y direction
    // const double Ly;
    // length in z direction
    // const double Lz;
    // number of cells in x direction
    // const int nx = 80;
    // number of cells in y direction
    // const int ny = 90;
    // number of cells in z direction
    // const int nz = 100;
    // x_vector
    // std::array<double,nx+1> x;
    // y_vector
    // std::array<double,ny+1> y;
    // z_vector
    // std::array<double,nz+1> z;
    // x_meshgrid
    // X;
    // y_meshgrid
    // Y;
    // z_meshgrid
    // Z;
    // [nx x ny x nz]-matrix with zeros inside and ones outside the boundary
    // boundary;
            
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
    // energy sharing in ionizing collision
    W;
    // maximum electron energy
    E_max;
    // maximal collision frequency:
    double nu_max;
    // collision counter
    unsigned int counter;
    // checks end of simulation: End =1 stops the simulation
    bool End = 0;
    // equilibrium time
    double T_sst = 0.0;
    // Counter of how many time steps there have been after steady state
    unsigned int count_sst = 0;

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
    std::array<MATRIX,PARTICLES_TYPES> r;
    // current velocity of electrons
    MATRIX v;
    // current acceleration of electrons
    MATRIX a;
    // current time-integrated velocity
    MATRIX v_int;
    // current time-integrated velocity-squared
    MATRIX v2_int;
    // collision indices for elastic collision
    ind_ela;
    // collision indices for excitation collision
    ind_exc;
    // collision indices for ionization collision
    ind_ion;
    // collision indices for attachment collision
    ind_att;
    // total number of all real collisions that happend
    int collisions = 0;
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
    std::vector<MeanData> mean; 
    // bulk transport data
    BulkData bulk;
    // flux transport data
    FluxData flux;
    // reaction rates
    RateDataConv rates_conv;
    RateDataCount rates_count;
    // energy data
    EnergyData E;
            
    // CLASS SOLVE_POISSON --> not needed here
    // charge density
    // rho;
    // electric potential
    // phi;
    // electric field function in x-direction
    double E_x;
    // electric field function in y-direction
    double E_y;
    // electric field function in z-direction
    double E_z;
            
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
    
    // Computes masses [kg] of the gas species
    void mass_in_kg();

    // Calculates the gas number density (in m^-3) by the ideal
    // gas law from pressure p and temperature Temp
    void gasNumberDensity();

    // Calculates absolute value of velocity and energy in eV
    std::pair<double,double> velocity2energy(const std::vector<double> & v);

    // Calculates maximal collision rate for gas mixture (in s^-1) 
    // void maximalCollFreq();                                // USA CLASSE "cross_sect"

    // Sets initial position and velocity of electrons
    void initialParticles();

    // Removes electrons outside the boundary 
    void surfaceInteraction();

    // Methods for the computation of the electric field:
    // particle2density(), solvePoisson_3D(), ...
    
    // Sets E parallel to z-direction --> replaces "solvePoissin_3D()"
    void set_E(); 

    // Generates random numbers p from an uniform distribution U[0,1];
    double random();
    std::vector<double> random(const unsigned N);

    // Performs non-collissional flight for electrons in electric field
    void freeFlight();

    // Gets electron collective data: mean position, mean
    // broadening in x,y and-direction, mean kinetic energy, electron number
    // and total electron current
    void collectMeanData();

    // Counts how many time steps there have been after steady state.
    void update_count_sst();

    // Calculates mean energy and EEDF data after steady state was reached
    void updateEnergyData();

    // Calculates flux data after steady state was reached
    void updateFluxData();

    // Calculates bulk data after steady state was reached
    void updateBulkData();

    // Calculates reaction rates after steady state was reached
    void updateReactionRates();
};

#endif 