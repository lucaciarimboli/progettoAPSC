#ifndef COLLISION_MATRIX_H
#define COLLISION_MATRIX_H

#include "CrossSectionsData.hpp"
#include <vector>
#include <random>
#include <algorithm>

class CollisionData {

public:
    // Constructor
    CollisionData(const int n_particles, const int n_react) :
        C(n_react, std::vector<double>(n_particles, 0.0)), Mass(n_react, 0.0), Loss(n_react, 0.0){}

    // Destructor
    ~CollisionData() = default;

    // Build collision matrix and compute indeces:
    void collisionmatrix(const int n_particles, const std::vector<double> & E_in_eV, const std::vector<double> & v_abs,
        const std::vector<double> & XS_energy, const std::vector<std::vector<table>> & XS,
        const std::vector<double> & mix, const std::vector<double> & mgas,
        const double E_max, const double nu_max, const double density){

        // Check if the number of particles matches the size of E_in_eV and v_abs:
        if( n_particles != E_in_eV.size() || n_particles != v_abs.size()){
            throw std::invalid_argument("Number of particles and corresponding velocities/energies do not match.");
        }

        // Initialize collision matrix with zeros:
        int n_react = C.size(); // Number of reactions
        C.resize(n_react, std::vector<double>(n_particles, 0.0));

        // Initialize mass and loss vectors with zeros:
        std::vector<double> mass(n_react, 0.0);
        std::vector<double> loss(n_react, 0.0);

        // Initialize collision indices:
        std::vector<size_t> col_ela;
        std::vector<size_t> col_exc;
        std::vector<size_t> col_ion;
        std::vector<size_t> col_att;

        // Note that the six vectors are needed just to store data
        // which will be used to fill the actual Loss, Mass and index vectors (that are class members).

        // Fill the collision matrix with the cross-section data:
        collisionfreq(mass, loss, E_in_eV, v_abs, XS_energy, XS, mix, mgas, E_max, nu_max, density,
            col_ela, col_exc, col_ion, col_att);

        // Simulate which interaction will happen for each particle:
        std::default_random_engine gen;
        std::uniform_real_distribution<> randu(0.0,1.0);
        std::vector<size_t> ind(n_particles,0);     // ind[i] will store the reaction index for i-esim electron
        for (size_t i = 0; i < n_particles; i++) {
            ind[i] = 1;
            double R = randu(gen); // [0,1] uniform distribution
            for( const std::vector<double> & c : C){
                if( R < c[i] ) ind[i]++;
            }
        }

        // Compute Mass and Loss vectors according to the collision indices:
        fill_Mass(ind, col_ela, mass);
        fill_Loss(ind, col_exc, loss);
        fill_Loss(ind, col_ion, loss);

        // Fill the vectors with collision indeces:
        find_collision_indeces(ind, mass, loss, col_ela, col_exc, col_ion, col_att);
    }


    // Getters:
    const double getC(int i, int j) const { return C[i][j]; }
    const std::vector<double>& getMass() const { return Mass; }
    const std::vector<double>& getLoss() const { return Loss; }
    const std::vector<size_t>& get_ind(std::string type) const {
        if (type == "ELASTIC") return ind_ela;
        else if (type == "EXCITATION") return ind_exc;
        else if (type == "IONIZATION") return ind_ion;
        else if (type == "ATTACHMENT") return ind_att;
        else throw std::invalid_argument("Invalid collision type");
    }

    // Count and return the number of real collisions happened:
    const int getCollisions() const {
        return ind_ela.size() + ind_exc.size() + ind_ion.size() + ind_att.size();
    }

private:
    std::vector<std::vector<double>> C;     // Collision matrix

    std::vector<double> Mass;               // Mass vector
    std::vector<double> Loss;               // Loss vector    
    
    std::vector<size_t> ind_ela;           // Collision indices for elastic collision
    std::vector<size_t> ind_exc;           // Collision indices for excitation collision
    std::vector<size_t> ind_ion;           // Collision indices for ionization collision
    std::vector<size_t> ind_att;           // Collision indices for attachment collision

    // Computes collision frequencies and fills "mass" and "loss" vectors:
    void collisionfreq( std::vector<double> & mass, std::vector<double> & loss,
        const std::vector<double> & E_in_eV, const std::vector<double> & abs_v,
        const std::vector<double> & XS_energy, const std::vector<std::vector<table>> & XS,
        const std::vector<double> & mix, const std::vector<double> & mgas,
        const double E_max, const double nu_max, const double density,
        std::vector<size_t> & col_ela, std::vector<size_t> & col_exc,
        std::vector<size_t> & col_ion, std::vector<size_t> & col_att) {
        
        // Copy energy grid to eventually manipulate it for interpolation without changing the original
        std::vector<double> x = XS_energy;
        size_t react_index = 0;

        bool flag_Emax = false;
        if( E_max > x.back()){
            flag_Emax = true;
            x.push_back(E_max);
        }
        
        for( size_t i = 0; i < XS.size(); i++){
            for( size_t j = 0; j < XS[i].size(); j++){
                if( XS[i][j].interact != EFFECTIVE){
                    // Same as for the energy grid, here for the cross-section data
                    std::vector<double> y = XS[i][j].section;
                
                    // Add data point at the energy E_max (if needed):
                    if(flag_Emax) y.push_back(y.back());

                    // Fill C by interpolation:
                    fill_C(x, y, E_in_eV, react_index, mix[i], density, abs_v, nu_max);

                    // Update the collision indices and Loss vector:
                    if( XS[i][j].interact == ELASTIC){
                        col_ela.push_back(react_index);
                        // loss[react_index] = 0.0;
                    } else if( XS[i][j].interact == EXCITATION){
                        col_exc.push_back(react_index);
                        loss[react_index] = XS[i][j].en_avg;        // VERIFICA SE E' CORRETTO !!
                    } else if( XS[i][j].interact == IONIZATION){
                        col_ion.push_back(react_index);
                        loss[react_index] = XS[i][j].en_avg;        // VERIFICA SE E' CORRETTO !!
                    } else if( XS[i][j].interact == ATTACHMENT){
                        col_att.push_back(react_index);
                        // loss[react_index] = 0.0;
                    }

                    // Update the mass vector:
                    mass[react_index] = mgas[i];
                    react_index++;
                }
            }
        }
    }

    // Fills the collision matrix C by interpolation and cumulative sum:
    void fill_C(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& xq,
        const size_t react_index, const double fract, double const density,
        const std::vector<double>& abs_v, double nu_max) {

        // This methods has been obtained from the "linear_interpolation" method used in other classes
        // by slighty modifying it to fill the collision matrix C instead of just performing the interpolation.

        if (x.empty() || y.empty() || x.size() != y.size()) {
            throw std::invalid_argument("Input vectors x and y must be non-empty and of the same size.");
        }

        for (size_t n = 0; n < xq.size(); n++) {
            double query = xq[n];
            size_t i = 0;

            if (query <= x.front()) {
            i = 0;
            } else if (query >= x.back()) {
            i = x.size() - 2;
            } else {
            while (i < x.size() - 1 && query > x[i + 1]) i++;
            }

            double t = (query - x[i]) / (x[i + 1] - x[i]);
            double sigma_f = y[i] + t * (y[i + 1] - y[i]);
            double c = fract * density * sigma_f * abs_v[n] / nu_max;

            // Update the collision matrix C by cumulative sum of the values obtained from the interpolated values
            if(react_index > 0) c += C[react_index-1][n];

            // Check if C is >1 or nan: if not update the collision matrix
            if(c > 1.0 || std::isnan(c))
                throw std::overflow_error("Collision matrix value exceeds 1 or is NaN");
            
            C[react_index][n] = c;
        }
    }

    void fill_Mass(const std::vector<size_t> ind, const std::vector<size_t>& coll_vect, const std::vector<double>& mass) {
        for( double col : coll_vect){
            auto it = std::find(ind.cbegin(), ind.cend(), col);
            // Check if the index is found:
            if (it == ind.cend())
                throw std::runtime_error("Collision index not found");
            // Fill the Mass/Loss vector at the found index:
            size_t index = std::distance(ind.cbegin(), it);
            Mass[index] = mass[col];
        }
    }

    void fill_Loss(const std::vector<size_t> & ind, const std::vector<size_t>& coll_vect, const std::vector<double>& loss) {
        for( double col : coll_vect){
            auto it = std::find(ind.cbegin(), ind.cend(), col);
            // Check if the index is found:
            if (it == ind.cend())
                throw std::runtime_error("Collision index not found");
            // Fill the Mass/Loss vector at the found index:
            size_t index = std::distance(ind.cbegin(), it);
            Loss[index] = loss[col];
        }
    }

    void find_collision_indeces( const std::vector<size_t>& ind,
        const std::vector<double>& mass, const std::vector<double>& loss,
        const std::vector<size_t>& col_ela, const std::vector<size_t>& col_exc,
        const std::vector<size_t>& col_ion, const std::vector<size_t>& col_att) {
        int col_ela_min = col_ela.front();
        int col_ela_max = col_ela.back();
        int col_exc_min = col_exc.front();
        int col_exc_max = col_exc.back();
        int col_ion_min = col_ion.front();
        int col_ion_max = col_ion.back();
        int col_att_min = col_att.front();
        int col_att_max = col_att.back();

        for (size_t i = 0; i < ind.size(); i++) {
            int val = ind[i];
            if (val >= col_ela_min && val <= col_ela_max) ind_ela.push_back(i);
            else if (val >= col_exc_min && val <= col_exc_max) ind_exc.push_back(i);
            else if (val >= col_ion_min && val <= col_ion_max) ind_ion.push_back(i);
            else if (val >= col_att_min && val <= col_att_max) ind_att.push_back(i);
        }
    }

};

#endif // COLLISION_MATRIX_Hx