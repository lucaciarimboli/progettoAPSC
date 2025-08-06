#include "utils/CollisionData.hpp"

CollisionData::CollisionData(const CrossSectionsData& Xsec, const std::vector<double>& mgas)
{
    const int n_react = Xsec.get_n_react(); // Number of reactions

    // Reserve space for vectors:
    mass.resize(n_react, 0.0);              
    loss.resize(n_react, 0.0);             
    col_ela.reserve(n_react);               
    col_exc.reserve(n_react);               
    col_ion.reserve(n_react);               
    col_att.reserve(n_react);    

    size_t react_index = 0;

    //for( size_t react_index = 0; react_index < n_react; react_index++){
    //    const table& tab = Xsec.get_full_xs_data()[react_index];
    for( const table& tab : Xsec.get_full_xs_data()){

        //if( tab.interact != mc::EFFECTIVE){

        // Update the collision indices and "loss" vector:
        if( tab.interact == mc::ELASTIC){
            col_ela.insert(react_index);
        } else if( tab.interact == mc::EXCITATION){
            col_exc.insert(react_index);
            loss[react_index] = tab.en_avg;        // VERIFICA SE E' CORRETTO !!
        } else if( tab.interact == mc::IONIZATION){
            col_ion.insert(react_index);
            loss[react_index] = tab.en_avg;        // VERIFICA SE E' CORRETTO !!
        } else if( tab.interact == mc::ATTACHMENT){
            col_att.insert(react_index);
        }

        // Update the "mass" vector:
        mass[react_index] = mgas[tab.specie_index];
        react_index++;
        //}
    }
}

const std::vector<size_t>& CollisionData::get_ind(mc::InteractionType type) const {
    if (type == mc::ELASTIC)         return ind_ela;
    else if (type == mc::EXCITATION) return ind_exc;
    else if (type == mc::IONIZATION) return ind_ion;
    else if (type == mc::ATTACHMENT) return ind_att;
    else throw std::invalid_argument("Invalid collision type");
}

const int CollisionData::getCollisions() const {
    // Count and return the number of real collisions happened:
    return ind_ela.size() + ind_exc.size() + ind_ion.size() + ind_att.size();
}

void CollisionData::ComputeIndeces(const int& n_electrons, const CrossSectionsData& Xsec,
    const std::vector<double>& E_in_eV, const std::vector<double>& v_abs,
    const std::vector<double>& mix, const double& density, const std::vector<double>& R)
{
    // Build collision matrix and compute indeces:

    // Define vector to store which reaction will happen for each electron:
    std::vector<size_t> ind(n_electrons);
    // ind.reserve(n_electrons);

    // Common factor in the collision matrix:
    const double common_factor = density / Xsec.get_nu_max();
    const std::vector<double>& XS_energy = Xsec.get_energy();
    const std::vector<table>& XS_data = Xsec.get_full_xs_data();

    std::vector<size_t> k(n_electrons);
    std::vector<double> t(n_electrons);
    #pragma omp parallel for schedule(static, 10000)
    for (int i = 0; i < n_electrons; i++) {
        auto it = std::upper_bound(XS_energy.cbegin(), XS_energy.cend(), E_in_eV[i]);
        const size_t idx = std::min(
            static_cast<size_t>(it - XS_energy.cbegin() - 1),
            XS_energy.size() - 2
        );
        k[i] = idx;
        t[i] = (E_in_eV[i] - XS_energy[idx]) / (XS_energy[idx + 1] - XS_energy[idx]);
    }

    // Compute the cumulative sum by electron of the collision matrix:
    #pragma omp parallel for schedule(static, 10000)
    for(int i = 0; i < n_electrons; i++){
        // Simulate which collision the i-esim electron undergoes:
        // ind.push_back(CollisionMatrix(R[i], XS_data, mix, t[i], k[i], v_abs[i]*common_factor));
        ind[i] = CollisionMatrix(R[i], XS_data, mix, t[i], k[i], v_abs[i]*common_factor);
    }

    // Compute Mass, Loss and indeces vectors:
    fill_Mass(ind);
    fill_Loss(ind);
    find_collision_indeces(ind);
}

const size_t CollisionData::CollisionMatrix(const double& R, const std::vector<table>& XS_data, const std::vector<double>& mix,
   const double& t, const size_t& k, const double& factor)
{    
    // Build collision matrix and update current index based on the random number
    double c = 0.0;             // cumsum of the collision matrix row corresponding to the current electron
    const size_t num_collisions = XS_data.size();

    for(size_t i = 0; i < num_collisions; i++){ 

        // Obtain by interpolation the xsec value corresponding to electron energy level:
        const std::vector<double>& yy = XS_data[i].section;
        const double sigma_f = yy[k] + t * (yy[k+1] - yy[k]);

        // Collision matrix element, corresponding to current electron, current interaction:
        c += factor * mix[XS_data[i].specie_index] * sigma_f;

        if( c > R) return i;
    }
    return num_collisions;  // this means that no collision has occourred.
}

void CollisionData::fill_Mass(const std::vector<size_t>& ind) {
    // Fill the Mass vector with the mass of the gas species

    // Initialize Mass vector with zeros:
    Mass.clear();
    Mass.resize(ind.size(), 0.0);

    for( size_t col : col_ela){
        auto it = std::find(ind.cbegin(), ind.cend(), col);
        // Check if the index is found:
        if (it != ind.cend()){
            // Fill the Mass vector at the found index:
            size_t index = std::distance(ind.cbegin(), it);
            Mass[index] = mass[col];
        }
    }
}

void CollisionData::fill_Loss(const std::vector<size_t> & ind) {
    // Fill the Loss vector with the energy loss of the gas species:

    // Initialize Loss vector with zeros:
    Loss.clear();
    Loss.resize(ind.size(), 0.0); 
    
  // Consider excitation collisions:
    for( size_t col : col_exc){
        auto it = std::find(ind.cbegin(), ind.cend(), col);
        // Check if the index is found:
        if (it != ind.cend()){
            // Fill the Loss vector at the found index:
            size_t index = std::distance(ind.cbegin(), it);
            Loss[index] = loss[col];
        }
    }

    // Consider ionization collisions:
    for( size_t col : col_ion){
        auto it = std::find(ind.cbegin(), ind.cend(), col);
        // Check if the index is found:
        if (it != ind.cend()){
        // Fill the Loss vector at the found index:
        size_t index = std::distance(ind.cbegin(), it);
        Loss[index] = loss[col];
        }
    }
}

void CollisionData::find_collision_indeces(const std::vector<size_t>& ind) {
    // Find the collision indeces for each type of collision:

    // Initialize index vectors as empty:
    ind_ela.clear(); ind_ela.reserve(ind.size());
    ind_exc.clear(); ind_exc.reserve(ind.size());
    ind_ion.clear(); ind_ion.reserve(ind.size());
    ind_att.clear(); ind_att.reserve(ind.size());

    for (size_t i = 0; i < ind.size(); i++) {
        const size_t val = ind[i];
        if (col_ela.count(val)) ind_ela.push_back(i);
        else if (col_exc.count(val)) ind_exc.push_back(i);
        else if (col_ion.count(val)) ind_ion.push_back(i);
        else if (col_att.count(val)) ind_att.push_back(i);
    }

    // Free not-used space in vectors:
    ind_ela.shrink_to_fit();
    ind_exc.shrink_to_fit();
    ind_ion.shrink_to_fit();
    ind_att.shrink_to_fit();
}
