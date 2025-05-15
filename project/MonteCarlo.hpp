#ifndef MONTECARLO_H
#define MONTECARLO_H

class CrossSectionsData;
class MeanData;
class EnergyData;
class FluxData;
class BulkData;
class RateDataCount;
class RateDataConv;
class CollisionData;

enum ParticleType {ELECTRONS = 0, CATIONS, ANIONS, PARTICLES_TYPES};
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

    // cross sections data
    CrossSectionsData Xsec;
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
    bool conserve = true;    
    // (1) isotropic, (0) non-isotropic scattering according to Vahedi et al.
    bool isotropic = true;
    // energy sharing in ionizing collision
    double W;
    // maximum electron energy
    double E_max;
    // maximal collision frequency
    double nu_max;
    // collision counter
    unsigned int counter;
    // checks end of simulation: End =1 stops the simulation
    bool End = false;
    // equilibrium time
    double T_sst = 0.0;
    // Counter of how many time steps there have been after steady state
    unsigned int count_sst = 0;

    // line number in output file:
    int line = 1;
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
    // Collision matrix and collision indeces:
    CollisionData C;
    // total number of all real collisions that happend
    int collisions = 0;
            
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
            
    // electric field function in x-direction
    double E_x;
    // electric field function in y-direction
    double E_y;
    // electric field function in z-direction
    double E_z;
            
    //Heat energy partition
    double EnergyLossElastic    = 0;
    double EnergyLossInelastic  = 0;
    double EnergyLossIonization = 0;
            
    //Status
    unsigned int converge = 0;

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
    std::pair<double,double> velocity2energy(const std::array<double,3> & v);

    // Calculates maximal collision rate for gas mixture (in s^-1) 
    void maximalCollFreq();

    // Sets initial position and velocity of electrons
    void initialParticles();

    // Removes electrons outside the boundary 
    void surfaceInteraction();
    
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

    // Decides which collision will happen for each electron
    void updateCollisionMatrix();

    // Perform collisions:
    void performCollisions();

    // Checks if the simulation has reached steady state
    void checkSteadyState();

    // Stops the simulation
    void endSimulation();

    private:

    // Computes cross product of two 3D vectors (needed in methods that perform collisions):
    std::array<double, 3> cross_product(const std::array<double, 3>& a, const std::array<double, 3>& b);
    // Performs elastic collision (isotropic or non-isotropic)
    void elasticCollision(const std::vector<size_t> & ind, const std::vector<double> & Mass);
    // Performs inelastic collision (isotropic or non-isotropic)
    void inelasticCollision(const std::vector<size_t> & ind, const std::vector<double> & Loss);
    // Performs ionization collision for electrons
    void ionizationCollision(const std::vector<size_t> & ind);
    // Performs attachment collision for electrons
    void attachmentCollision(const std::vector<size_t> & ind);

};

#endif 