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

        if( tab.interact != mc::EFFECTIVE){

            // Update the collision indices and "loss" vector:
            if( tab.interact == mc::ELASTIC){
                col_ela.insert(react_index);
                // loss[react_index] = 0.0;
            } else if( tab.interact == mc::EXCITATION){
                col_exc.insert(react_index);
                loss[react_index] = tab.en_avg;        // VERIFICA SE E' CORRETTO !!
            } else if( tab.interact == mc::IONIZATION){
                col_ion.insert(react_index);
                loss[react_index] = tab.en_avg;        // VERIFICA SE E' CORRETTO !!
            } else if( tab.interact == mc::ATTACHMENT){
                col_att.insert(react_index);
               // loss[react_index] = 0.0;
            }

            // Update the "mass" vector:
            mass[react_index] = mgas[tab.specie_index];
            react_index++;
        }
    }
}

void CollisionData::ComputeIndeces(const int n_electrons, const CrossSectionsData& Xsec,
    const std::vector<double>& E_in_eV, const std::vector<double>& mix, const double density, const std::vector<double>& R)
{
    // Check if the number of particles matches the size of E_in_eV and v_abs:
    //if( n_electrons != E_in_eV.size()){
    //    throw std::invalid_argument("Number of particles and corresponding velocities/energies do not match.");
    //}

    // Define vector to store which reaction will happen for each electron:
    std::vector<size_t> ind;
    ind.resize(n_electrons,1);

    // Compute the cumulative sum by electron of the collision matrix:
    for( int e_index = 0; e_index < n_electrons; e_index++){

        // Compute the collision matrix row corresponding to the current electron:
        //const size_t collision_index = CollisionMatrix(R[e_index], Xsec, mix, E_in_eV[e_index], density);
        const size_t collision_index = CollisionMatrix(R[e_index], Xsec, mix, E_in_eV[e_index], density);

        // Assign the collision index to the corresponding electron:
        ind[e_index] = collision_index;
    }

    // Compute Mass, Loss and indeces vectors:
    fill_Mass(ind);
    fill_Loss(ind);
    find_collision_indeces(ind);
}

const size_t CollisionData::CollisionMatrix(const double R, const CrossSectionsData& Xsec,
    const std::vector<double>& mix, const double E_in_eV, const double density)
{    
    const std::vector<double>& xx = Xsec.get_energy();  // energy bins (eV) from xs data
    double nu_max = Xsec.get_nu_max();

    // Compute the absolute value of the velocity:
    double v_abs = std::sqrt(2.0 * E_in_eV * mc::q0 / mc::me); // m/s

    size_t react_index = 0;     // Index for the current reaction
    size_t collision_index = 0; // Index for the xollision which the electron is undergoing

    double c = 0.0;             // cumsum of the collision matrix row corresponding to the current electron

    for(const table& tab : Xsec.get_full_xs_data()){ 

        if( tab.interact != mc::EFFECTIVE){

            const std::vector<double>& yy = tab.section; // cross-section data

            // Compute the collision frequency by interpolation:
            size_t k = 0;
            if (E_in_eV <= xx.front()) {
                k = 0;
            } else if (E_in_eV >= xx.back()) {
                k = xx.size() - 2;
            } else {
                while (k < xx.size() - 1 && E_in_eV > xx[k + 1]) k++;
            }
            double t = (E_in_eV - xx[k]) / (xx[k + 1] - xx[k]);
            double sigma_f = yy[k] + t * (yy[k + 1] - yy[k]);
            c += mix[tab.specie_index] * density * sigma_f * v_abs / nu_max;

            // Check if C is >1 or nan: if not update the collision matrix
            /*if(std::isnan(c)){
                throw std::overflow_error("Collision matrix value is NaN"); 
            }*/
            
            // THE COMMENTED CODE BELOW WORKS BUT TOTAL NUMBER OF COLLISIONS EXCEEDS NUMBER OF ELECTRONS!
            /*if(c > 1.0) break;

            // Electron undergoes the first reaction "i" s.t. C[i] > R:
            if( R > c) collision_index++;

            react_index++;*/


            // NEW (MAYBE MORE EFFCIENT) VERSION:
            if( c > R) break;
            collision_index++;
            react_index++;
        }
    }
    //return (collision_index == 0) ? 0 : (collision_index - 1);
    return collision_index;
}

void CollisionData::fill_Mass(const std::vector<size_t>& ind) {

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
        }// else {
        //    throw std::runtime_error("Collision index not found");
        //}
    }

    // In the data provided, there are no elastic collisions!!!
    // This means that the Mass vector would be all zeros and the ind_ela must be empty!
}

void CollisionData::fill_Loss(const std::vector<size_t> & ind) {

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
        }// else {
        //    throw std::runtime_error("Collision index not found");
        //}
    }

    // Consider ionization collisions:
    for( size_t col : col_ion){
        auto it = std::find(ind.cbegin(), ind.cend(), col);
        // Check if the index is found:
        if (it != ind.cend()){
        // Fill the Loss vector at the found index:
        size_t index = std::distance(ind.cbegin(), it);
        Loss[index] = loss[col];
        }// else {
        //    throw std::runtime_error("Collision index not found");
        //}
    }
}

void CollisionData::find_collision_indeces(const std::vector<size_t>& ind) {
    // Initialize index vectors as empty:
    ind_ela.clear(); ind_ela.reserve(ind.size());
    ind_exc.clear(); ind_exc.reserve(ind.size());
    ind_ion.clear(); ind_ion.reserve(ind.size());
    ind_att.clear(); ind_att.reserve(ind.size());

    for (size_t i = 0; i < ind.size(); i++) {
        size_t val = ind[i];
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