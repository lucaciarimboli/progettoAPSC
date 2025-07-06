#include "core/MonteCarlo.hpp"

MonteCarlo::MonteCarlo( const std::vector<std::string> & gas, const std::vector<double> & mix,
                const double EN, const double p, const double T, const unsigned N0, const unsigned Ne_max,
                const double W, const double E_max, const double dE,
                const double w_err, const double DN_err, const unsigned col_equ, const unsigned col_max,
                const std::array<double,3> & pos_xyz, const std::array<double,3> & sigma_xyz,
                const bool conserve, const bool isotropic):

    N0(N0), N(p/(mc::kB * T)), gas(gas), mix(mix), Xsec(gas, E_max, mix, N),
    w_err(std::abs(w_err)), DN_err(std::abs(DN_err)), Ne_max(Ne_max), col_equ(col_equ), col_max(col_max),
    conserve(conserve), isotropic(isotropic), W(W), E_max(E_max), EN(EN), t(1,0.0), dt(0.0), v(N0, {0.0, 0.0, 0.0}),
    v_int(N0, {0.0, 0.0, 0.0}), v2_int(N0, {0.0, 0.0, 0.0}), mean(1, MeanData(pos_xyz, sigma_xyz, N0)),
    E([&]() {
        std::vector<double> energy_bins;
        energy_bins.reserve(static_cast<size_t>(E_max / dE) + 1);
        for (double energy = 0.0; energy <= E_max; energy += dE) 
            energy_bins.push_back(energy);
        return EnergyData(energy_bins);
    }()),
    bulk(), flux(), rates_conv(Xsec, E, mix), rates_count(N, conserve),
    gen(std::random_device{}()), randu(0.0,1.0) {
        
    // The check for the validity of the gas species is done in "CrossSectionsData" constructor

    // Check mix vector validity:
    checkFractionSum();

    // Compute mass of the gas species in kg:
    mass_in_kg();

    // Initialize energy data:
    //std::vector<double> energy_bins;
    //energy_bins.reserve(static_cast<size_t>(E_max / dE) + 1);  // Avoid reallocations
    //for (double E = 0.0; E <= E_max; E += dE) energy_bins.push_back(E);
    //E = EnergyData(energy_bins);
    // E = EnergyData(E_max, dE);

    // Set electric field E (constant and uniform):
    const double E_x = 0;
    const double E_y = 0;
    const double E_z = EN * N * 1e-21; // Electric field strength in V/m, EN is the energy in eV

    // Initialize acceleration array after the electric field was set:
    // (since E is constant and uniform, a does not change and is the same for every electron)
    a = {mc::q0/mc::me*E_x, mc::q0/mc::me*E_y, mc::q0/mc::me*E_z};

    // Initialize particles by Gaussian distribution:
    initialParticles(pos_xyz, sigma_xyz);

    // Initialize data for computation of collision frequencies:
    C = CollisionData(Xsec, mgas);
}

void MonteCarlo::checkFractionSum(){
    // Checks if sum of gas fractions is equal to 1
    // if not the case: last entry of mix will be corrected

    // Check if the size of gas and mix vectors are equal
    if(gas.size() != mix.size()){
        std::cerr << "Error: gas and mix vectors have different sizes!" << std::endl;
        return;
    }

    // Check if the sum of mixing ratios is equal to 1
    // If not, correct the last entry of mix
    double sumMix = std::accumulate(mix.cbegin(),mix.cend(),0.0);        
    if(sumMix != 1){
            mix.back() = 1 - std::accumulate(mix.cbegin(),mix.cend()-1,0.0);
            std::cout << "Sum of mixing ratios (in mix) NOT equal to one!" << std::endl;
    }
}

void MonteCarlo::mass_in_kg(){
    // Computes mass of the gas species in kg
    mgas.reserve(gas.size());

    MolMass mol_mass;
    for(const std::string & substance : gas){
        mol_mass.set_substance(substance);
        mol_mass.Compute_M();
        const double M = mol_mass.get_front_M();
        mgas.push_back( M / (mc::Na * 1000) );
    }
}

std::pair<double,double> MonteCarlo::velocity2energy(const std::array<double,3> & v) const {
    // calculates absolute value of velocity abs_v and energy
    // E_in_eV in eV for one electron.
    double abs_v = std::sqrt(std::inner_product(v.cbegin(), v.cend(), v.cbegin(), 0.0));
    double E_in_eV =  0.5 * mc::me *  abs_v * abs_v / mc::q0;
    return std::make_pair(abs_v,E_in_eV);
}

void MonteCarlo::initialParticles(const std::array<double,3> & pos_xyz, const std::array<double,3> & sigma_xyz){
    // Sets initial position and velocity of electrons

    // Allocate memory to electron population positions:
    r[mc::ELECTRONS].clear();
    r[mc::ELECTRONS].reserve(N0);

    // Initialize cations and anions populations as empty:
    r[mc::CATIONS].clear();
    r[mc::ANIONS].clear();

    // Standard Gaussian distribution
    std::normal_distribution<double> randn(0.0,1.0);

    for (size_t n = 0; n < N0; n++) {
        // Initialize to pos_xyz + noise if sigma_xyz != 0
        r[mc::ELECTRONS].emplace_back(std::array<double, 3>{
            pos_xyz[0] + sigma_xyz[0] * randn(gen),
            pos_xyz[1] + sigma_xyz[1] * randn(gen),
            pos_xyz[2] + sigma_xyz[2] * randn(gen)
        });
    }

    // Electrons velocities are initialized with 0.0, add noise:
    // for( auto &vel : v){
    //     for( size_t i=0; i<=2; i++){
    //         vel[i] += 1e-6 * randu(gen);
    //     }
    // }
}

void MonteCarlo::freeFlight(){
    // Performs non-collisional flight for electrons in electric field
    dt = - std::log(randu(gen)) / Xsec.get_nu_max(); // generates time step
    const double dt2 = dt * dt;
    const double dt3 = dt2 * dt;

    // Update vector time:
    t.push_back(t.back() + dt);
    
    // Number of electrons:
    const double ne = v.size();

    if(T_sst > 0.0){
        count_sst++;
        v_int.clear();
        v2_int.clear();
        v_int.reserve(ne);
        v2_int.reserve(ne);
        for(size_t i = 0; i < ne; i++){
            // integrated velocity
            v_int.emplace_back(std::array<double, 3>{
                v[i][0]*dt + 0.5 * a[0] * dt2,
                v[i][1]*dt + 0.5 * a[1] * dt2,
                v[i][2]*dt + 0.5 * a[2] * dt2
            });

            // integrated velocity squared
            v2_int.emplace_back(std::array<double, 3>{
                v[i][0]*v[i][0]*dt + a[0]*v[i][0] * dt2 + a[0]*a[0] / 3 * dt3,
                v[i][1]*v[i][1]*dt + a[1]*v[i][1] * dt2 + a[1]*a[1] / 3 * dt3,
                v[i][2]*v[i][2]*dt + a[2]*v[i][2] * dt2 + a[2]*a[2] / 3 * dt3
            });
        }
    }

    // Update space and velocity:
    for(size_t i = 0; i < ne; i++){
        r[mc::ELECTRONS][i][0] += v[i][0]*dt + 0.5 * a[0] * dt2;
        v[i][0] += a[0] * dt;

        r[mc::ELECTRONS][i][1] += v[i][1]*dt + 0.5 * a[1] * dt2;
        v[i][1] += a[1] * dt;

        r[mc::ELECTRONS][i][2] += v[i][2]*dt + 0.5 * a[2] * dt2;
        v[i][2] += a[2] * dt;
    }


    /*
    // Con "-" invece di "+" perché l'elettrone si muove in direzione opposta al campo E!
    for(size_t i = 0; i < ne; i++){
        r[mc::ELECTRONS][i][0] += v[i][0]*dt - 0.5 * a[0] * dt2;
        v[i][0] -= a[0] * dt;

        r[mc::ELECTRONS][i][1] += v[i][1]*dt - 0.5 * a[1] * dt2;
        v[i][1] -= a[1] * dt;

        r[mc::ELECTRONS][i][2] += v[i][2]*dt - 0.5 * a[2] * dt2;
        v[i][2] -= a[2] * dt;
    }
    */

}

void MonteCarlo::collectMeanData(){
    // Gets electron collective data: mean position, mean
    // broadening in x,y and-direction, mean kinetic energy, electron number
    // and total electron current

    // Update mean vector for the new time step
    mean.emplace_back(mean.back().get_particles(), r[mc::ELECTRONS], v);
}

void MonteCarlo::updateEnergyData(){
    // Calculates mean energy and EEDF data after steady state was reached

    const unsigned int ne = v.size(); // number of electrons
    t_total += dt * ne; // sum of all times for all electrons
    E.mean_energy(v2_int, t_total);

    // Compute kinetic energy for each electron:
    std::vector<double> E_in_eV(ne);    
    std::transform(v.begin(), v.end(), E_in_eV.begin(), [this](const std::array<double, 3>& vi) {
        return velocity2energy(vi).second;
    });
    E.energy_bins(E_in_eV);
    E.compute_distribution_function();  
}

void MonteCarlo::updateFluxData(){
    // Calculates flux data after steady state was reached
    flux.compute_drift_velocity(v_int,t_total);
    flux.compute_diffusion_const(r[mc::ELECTRONS],v,N, count_sst);
}

void MonteCarlo::updateBulkData(){
    // Calculates bulk data after steady state was reached
    bulk.update_bulk(t,count_sst,mean,N);
}

void MonteCarlo::updateReactionRates(){
    // Calculates reaction rates after steady state was reached

    // 1. REACTION RATES BY COUNTING:
    rates_count.setTime(t,count_sst);
    rates_count.setParticles(mean,count_sst);
    rates_count.computeRates();

    // 2. REACTION RATES BY CONVOLUTION:
    rates_conv.computeRates();
};

void MonteCarlo::updateCollisionMatrix(){
    // Decides which collision will happen for each electron

    // Extract number of electrons:
    const int num_particles = v.size();

    // From electron velocity compute energies in eV and extract the velocity modules:
    std::vector<double> E_in_eV;
    //std::vector<double> v_abs;
    E_in_eV.reserve(num_particles);
    //v_abs.reserve(num_particles);
    std::transform(v.begin(), v.end(), std::back_inserter(E_in_eV), [this](const std::array<double, 3>& vi) {
            return velocity2energy(vi).second;
        });
    /*std::transform(v.begin(), v.end(), std::back_inserter(v_abs), [this](const std::array<double, 3>& vi) {
            return velocity2energy(vi).first;
        });*/

    // Generate random numbers for collision indeces:
    std::vector<double> R(num_particles);
    std::generate(R.begin(), R.end(), [this]() { return randu(gen); });
    
    // Build collision matrix and compute indeces:
    C.ComputeIndeces(num_particles, Xsec, E_in_eV, mix, N, R); // I express v_abs in terms of E_in_eV inside the computations to avoid allocating memory for it
    // Update total number of collisions:
    collisions += C.getCollisions();
}

void MonteCarlo::performCollisions(){
    // Perform collisions according to the type of collision

    const std::vector<double> & Mass = C.getMass();
    const std::vector<double> & Loss = C.getLoss();

    const std::vector<size_t> & ind_ela = C.get_ind(mc::ELASTIC);
    const std::vector<size_t> & ind_exc = C.get_ind(mc::EXCITATION);
    const std::vector<size_t> & ind_ion = C.get_ind(mc::IONIZATION);
    const std::vector<size_t> & ind_att = C.get_ind(mc::ATTACHMENT);

    // If any, perform elastic collisions:
    if(!(ind_ela.empty())){
        elasticCollision(ind_ela, Mass);
    }
    // If any, perform inelastic collisions:
    if(!(ind_exc.empty())){
        inelasticCollision(ind_exc, Loss);
    }
    // if any, perform ionization collisions:
    if(!(ind_ion.empty())){
        ionizationCollision(ind_ion, Loss);
    }
    // If any, perform attachment collision:
    if(!(ind_att.empty())){
        attachmentCollision(ind_att);
    }
}

std::array<double, 3> MonteCarlo::cross_product(const std::array<double, 3>& a, const std::array<double, 3>& b) const {
    // Function to compute the cross product of two 3D vectors
    return {
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
}

void MonteCarlo::elasticCollision(const std::vector<size_t> & ind, const std::vector<double> & Mass){
    // Performs elastic collision (isotropic or non-isotropic)

    double E_1 = 0.0;                       // total energy before collision
    double E_2 = 0.0;                       // total energy after collision
    std::array<double,3> e_x = {1.0, 0.0, 0.0};  // x-direction versor

    double sin_phi;
    double cos_phi;
    double sin_xsi;
    double cos_xsi;
    double sin_theta;
    double cos_theta;

    for(size_t i = 0; i < ind.size(); i++){

        size_t el_index = ind[i];       // electron index
        std::pair<double,double> v2e = velocity2energy(v[el_index]);
        // Compute energy before collision
        E_1 += v2e.second;

        // Compute incident direction of scattering electrons:
        std::array<double,3> e_1 = {v[el_index][0]/v2e.first, v[el_index][1]/v2e.first, v[el_index][2]/v2e.first};

        // Randomly generate phi: azimuthal angle
        double phi = 2 * M_PI * randu(gen);
        sin_phi = std::sin(phi);
        cos_phi = std::cos(phi);

        // Randomly generate xsi: electron scattering angle
        if(isotropic) cos_xsi = 1 - 2 * randu(gen);
        else cos_xsi = (2 + v2e.second - 2 * (1+v2e.second) * randu(gen)) / v2e.second;
        sin_xsi = std::sqrt(1 - cos_xsi * cos_xsi);

        // Compute theta: angle between x-axis and incident velocity:
        cos_theta = e_1[0];
        sin_theta = std::sqrt(1 - cos_theta * cos_theta);

        // Compute the new direction e_2 of the scattered electron:
        std::array<double,3> cross1 = cross_product(e_1, e_x);
        std::array<double,3> cross2 = cross_product(e_x, e_1);
        std::array<double,3> cross3 = cross_product(e_1, cross2);
        std::array<double, 3> e_2 = {0.0, 0.0, 0.0};

        // ( Avoid division by zero in case theta is very small):
        double inverse_sin_theta = (sin_theta > 1e-10) ? 1.0 / sin_theta : 0.0;
        for (int j = 0; j < 3; j++) {
            e_2[j] = cos_xsi * e_1[j] +
                    sin_xsi * sin_phi * inverse_sin_theta * cross1[j] +
                    sin_xsi * cos_phi * inverse_sin_theta * cross3[j];
        }

        // Normalize e_2:
        double norm = std::sqrt(e_2[0]*e_2[0] + e_2[1]*e_2[1] + e_2[2]*e_2[2]);
        for (int j = 0; j < 3; ++j) e_2[j] /= norm;

        // Compute energy after the elastic collision:
        double E_2_el = std::max(0.0, v2e.second*(1 - 2*mc::me/Mass[el_index] * (1 - cos_xsi)));
        E_2 += E_2_el;

        // Update velocity:
        double v2_abs = std::sqrt(2.0 * E_2_el * mc::q0 / mc::me);
        for (int j = 0; j < 3; j++) {
            v[el_index][j] = v2_abs * e_2[j];
        }
    }

    // Update total energy loss counter:
    EnergyLossElastic += E_1 - E_2;
}

void MonteCarlo::inelasticCollision(const std::vector<size_t> & ind, const std::vector<double> & Loss){
    // Performs inelastic collision (isotropic or non-isotropic)

    double E_1 = 0.0;                       // total energy before collision
    double E_2 = 0.0;                       // total energy after collision
    std::array<double,3> e_x = {1.0, 0.0, 0.0};  // x-direction versor

    double sin_phi;
    double cos_phi;
    double sin_xsi;
    double cos_xsi;
    double sin_theta;
    double cos_theta;

    for(size_t i = 0; i < ind.size(); i++){

        size_t el_index = ind[i];       // electron index
        std::pair<double,double> v2e = velocity2energy(v[el_index]);
        // Compute energy before collision
        E_1 += v2e.second;

        // Compute incident direction of scattering electrons:
        std::array<double,3> e_1 = {v[el_index][0]/v2e.first, v[el_index][1]/v2e.first, v[el_index][2]/v2e.first};

        // Randomly generate phi: azimuthal angle
        double phi = 2 * M_PI * randu(gen);
        sin_phi = std::sin(phi);
        cos_phi = std::cos(phi);

        // Randomly generate xsi: electron scattering angle
        if(isotropic) cos_xsi = 1 - 2 * randu(gen);
        else cos_xsi = (2 + v2e.second - 2 * (1+v2e.second) * randu(gen)) / v2e.second;
        sin_xsi = std::sqrt(1 - cos_xsi * cos_xsi);

        // Compute theta: angle between x-axis and incident velocity:
        cos_theta = e_1[0];
        sin_theta = std::sqrt(1 - cos_theta * cos_theta);

        // Compute the new direction e_2 of the scattered electron:
        std::array<double,3> cross1 = cross_product(e_1, e_x);
        std::array<double,3> cross2 = cross_product(e_x, e_1);
        std::array<double,3> cross3 = cross_product(e_1, cross2);
        std::array<double, 3> e_2 = {0.0, 0.0, 0.0};

        // ( Avoid division by zero in case theta is very small):
        double inverse_sin_theta = (sin_theta > 1e-10) ? 1.0 / sin_theta : 0.0;
        for (int j = 0; j < 3; j++) {
            e_2[j] = cos_xsi * e_1[j] +
                sin_xsi * sin_phi * inverse_sin_theta * cross1[j] +
                sin_xsi * cos_phi * inverse_sin_theta * cross3[j];
        }

        // Normalize e_2:
        double norm = std::sqrt(e_2[0]*e_2[0] + e_2[1]*e_2[1] + e_2[2]*e_2[2]);
        for (int j = 0; j < 3; ++j) e_2[j] /= norm;

        // Compute energy after the elastic collision:
        double E_2_el = std::max(0.0, v2e.second - Loss[el_index]);
        E_2 += E_2_el;

        // Update velocity:
        double v2_abs = std::sqrt(2.0 * E_2_el * mc::q0 / mc::me);
        for (int j = 0; j < 3; j++) {
            v[el_index][j] = v2_abs * e_2[j];
        }
    }    

    // Update total energy loss counter:
    EnergyLossInelastic += E_1 - E_2;
}

void MonteCarlo::ionizationCollision(const std::vector<size_t> & ind, const std::vector<double> & Loss){
    // Performs ionization collision for electrons

    double E_1 = 0.0;                       // total energy before collision
    double E_2 = 0.0;                       // total energy after collision
    std::array<double,3> e_x = {1.0, 0.0, 0.0};  // x-direction versor

    // Number of created electrons by ionization:
    int delta_Ne = ind.size();

    double sin_phi;
    double cos_phi;
    double sin_xsi;
    double cos_xsi;
    double sin_theta;
    double cos_theta;

    for(int i = 0; i < delta_Ne; i++){

        size_t el_index = ind[i];       // electron index
        std::pair<double,double> v2e = velocity2energy(v[el_index]);
        // Compute energy before collision
        E_1 += v2e.second;

        // Compute incident direction of scattering electrons:
        std::array<double,3> e_1 = {v[el_index][0]/v2e.first, v[el_index][1]/v2e.first, v[el_index][2]/v2e.first};

        // Randomly generate phi: azimuthal angle
        double phi = 2 * M_PI * randu(gen);
        sin_phi = std::sin(phi);
        cos_phi = std::cos(phi);

        // Randomly generate xsi: electron scattering angle
        if(isotropic) cos_xsi = 1 - 2 * randu(gen);
        else cos_xsi = (2 + v2e.second - 2 * (1+v2e.second) * randu(gen)) / v2e.second;
        sin_xsi = std::sqrt(1 - cos_xsi * cos_xsi);

        // Compute theta: angle between x-axis and incident velocity:
        cos_theta = e_1[0];
        sin_theta = std::sqrt(1 - cos_theta * cos_theta);

        // Compute the new direction e_2 of the scattered electron:
        std::array<double,3> cross1 = cross_product(e_1, e_x);
        std::array<double,3> cross2 = cross_product(e_x, e_1);
        std::array<double,3> cross3 = cross_product(e_1, cross2);
        std::array<double, 3> e_2 = {0.0, 0.0, 0.0};

        // ( Avoid division by zero in case theta is very small):
        double inverse_sin_theta = (sin_theta > 1e-10) ? 1.0 / sin_theta : 0.0;
        for (int j = 0; j < 3; j++) {
            e_2[j] = cos_xsi * e_1[j] +
                    sin_xsi * sin_phi * inverse_sin_theta * cross1[j] +
                    sin_xsi * cos_phi * inverse_sin_theta * cross3[j];
        }

        // Normalize e_2:
        double norm = std::sqrt(e_2[0]*e_2[0] + e_2[1]*e_2[1] + e_2[2]*e_2[2]);
        for (int j = 0; j < 3; ++j) e_2[j] /= norm;

        // Compute energy after the elastic collision:
        double E_2_el = std::max(0.0, v2e.second - Loss[el_index]);
        E_2 += E_2_el;

        // Update velocity:
        double v2_abs = std::sqrt(2.0 * W * E_2_el * mc::q0 / mc::me); // W is the energy sharing in ionizing collision
        for (int j = 0; j < 3; j++) {
            v[el_index][j] = v2_abs * e_2[j];
        }

        // Add the new electron to the simulation:
        double v_abs_newe = std::sqrt(2.0 * (1-W) * E_2 * mc::q0 / mc::me);
        v.push_back({
            - v_abs_newe * e_2[0],
            - v_abs_newe * e_2[1],
            - v_abs_newe * e_2[2]
        });
        r[mc::ELECTRONS].push_back({
            r[mc::ELECTRONS][el_index][0],
            r[mc::ELECTRONS][el_index][1],
            r[mc::ELECTRONS][el_index][2]
        });

        // Add the new cation:
        r[mc::CATIONS].push_back({
            r[mc::ELECTRONS][el_index][0],
            r[mc::ELECTRONS][el_index][1],
            r[mc::ELECTRONS][el_index][2]
        });
    }

    // If required, enforce electron population conservation:
    if(conserve){
        for( int i = 0; i < delta_Ne; i++){
            // select a random electron:
            int random_index = static_cast<int>(randu(gen) * (r[mc::ELECTRONS].size()));
            // remove it from the simulation:
            r[mc::ELECTRONS].erase(r[mc::ELECTRONS].begin() + random_index);
            v.erase(v.begin() + random_index);
        }
    }
    
    // Update number of particles in mean data:
    mean.back().add_new_particles({delta_Ne, delta_Ne, 0});

    // Update total energy loss counter:
    EnergyLossIonization += E_1 - E_2;
}

void MonteCarlo::attachmentCollision(const std::vector<size_t> & ind){
    // Performs attachment collision for electrons

    // Number of created anions / removed electrons by attachment:
    int delta_Ne = ind.size();

    // Note: "ind" is built to be sorted in ascending order (see CollsionData::find_collision_indeces).
    // To remove safely the electrons without changing the other indeces, the loop is performed in reverse order.
    for (auto it = ind.rbegin(); it != ind.rend(); it++) {
    
        size_t el_index = *it;       // electron index

        // Add the new anion to the simulation:
        r[mc::ANIONS].push_back({
            r[mc::ELECTRONS][el_index][0],
            r[mc::ELECTRONS][el_index][1],
            r[mc::ELECTRONS][el_index][2]
        });

        // Remove the electron from the simulation:
        r[mc::ELECTRONS].erase(r[mc::ELECTRONS].begin() + el_index);
        v.erase(v.begin() + el_index);

        // If required, enforce electron population conservation:
        if(conserve){
            // select a random electron:
            int random_index = static_cast<int>(randu(gen) * (r[mc::ELECTRONS].size()));
            // clone it to compensate the removed one:
            r[mc::ELECTRONS].push_back(r[mc::ELECTRONS][random_index]);
            v.push_back(v[random_index]);
        }
    }

    // Update number of particles in mean data:
    mean.back().add_new_particles({- delta_Ne, 0, delta_Ne});
}

void MonteCarlo::checkSteadyState(){
    // Checks if the simulation has reached steady state.
    // If so, it updates the time step and the number of collisions

    // If sst has not been reached yet, check after last iteration:
    if(count_sst == 0 && collisions/1e6 >= line && collisions >= col_equ){

        // check if the interval 80-90 % of energy data is larger than the interval 90-100%:
        int N = mean.size();
        int n = std::round( N / 10.0);

        double sum1 = 0.0, sum2 = 0.0;
        for (int i = N - 2 * n; i < N - n; i++) {
            sum1 += mean[i].get_energy();
        }
        for (int i = N - n; i < N; i++) {
            sum2 += mean[i].get_energy();
        }

        if(sum1 >= sum2){
            // Steady state has been reached
            T_sst = t.back();
            std::cout << "\n Steady state reached at t = " << T_sst * 1e9 << " ns\n";
            std::cout << " Number of iterations until steady state: " << t.size()-1 << "\n";
            std::cout << " Number of collisions until steady state: " << collisions << "\n";

            // For debugging purposes:
            std::cout << " Size of r[ELECTRONS] at steady state: " << v.size() << "\n";

            collisions = 0;
            line = 1;
        }
    }
}

bool MonteCarlo::endSimulation() {
    // Stops the simulation

    // End simulation if too many electrons
    if (!conserve && v.size() > Ne_max) {
        converge = 1;
        std::cout << "\n Simulation ended: maximum number of electrons reached\n";
        return true;
    }

    // End simulation if no electrons
    if (!conserve && v.empty()) {
        converge = 2;
        std::cout << "\n Simulation ended: number of electrons is zero\n";
        return true;
    }

    // End simulation if relative errors in w and D are below thresholds
    if (!(bulk.is_empty())) {
        const std::array<double, 3> & w_bulk = bulk.get_w();
        const std::array<double, 3> & w_bulk_err = bulk.get_w_err();
        const std::array<double, 3> & DN_bulk = bulk.get_DN();
        const std::array<double, 3> & DN_bulk_err = bulk.get_DN_err();  
        if (std::abs(w_bulk_err[2] / w_bulk[2]) < w_err) {
            //&& std::abs(DN_bulk_err[2] / DN_bulk[2]) < DN_err) {

            converge = 0;
            std::cout << "\n Simulation ended: errors in w < " 
                      << 100 * w_err << "% and D < "
                      << 100 * DN_err << "%\n";
            return true;
        }
    }

    // End simulation if number of collisions exceeds maximum
    if (collisions > col_max) {
        converge = 3;
        std::cout << "\n Simulation ended: maximum number of collisions reached\n";
        return true;
    }

    return false;
}

void MonteCarlo::printOnScreen() {
    // Prints the status of relevant variables every 10^6 collisions

    if ((collisions / 1e6) >= line) {
        line += 1;

        if (!(bulk.is_empty())) {

            const std::array<double, 3> & w_bulk = bulk.get_w();
            const std::array<double, 3> & w_bulk_err = bulk.get_w_err();
            const std::array<double, 3> & DN_bulk = bulk.get_DN();
            const std::array<double, 3> & DN_bulk_err = bulk.get_DN_err();  
            const std::array<double, 3> & w_flux = flux.get_w();
            const std::array<double, 3> & DN_flux = flux.get_DN();
 
            std::printf(
                "\n Relative Error Bulk Velocity: %.3e\n"
                " Relative Error Diffusion Coeff: %.3e\n"
                " Collisions: %lu\n"
                " Electrons: %i\n"
                " Mean Energy: %.3e eV\n"
                " Bulk Velocity: %.3e m/s\n"
                " Flux Velocity: %.3e m/s\n"
                " Bulk Diffusion Coeff: %.2e (ms)^-1\n"
                " Flux Diffusion Coeff: %.2e (ms)^-1\n"
                " Effective React Rate (counted): %.2e m^3/s\n"
                " Effective React Rate (computed): %.2e m^3/s\n"
                " Alpha: %.3e m^-1\n"
                " Eta: %.3e m^-1\n",

                std::abs(w_bulk_err[2] / w_bulk[2]),
                std::abs(DN_bulk_err[2] / DN_bulk[2]),
                collisions,
                mean.back().get_particles()[mc::ELECTRONS],
                E.get_E_mean(),
                w_bulk[2],
                w_flux[2],
                DN_bulk[2],
                DN_flux[2],
                rates_count.getRate(mc::EFFECTIVE),
                rates_conv.getRate(mc::EFFECTIVE),
                rates_count.getRate(mc::IONIZATION) * 2.4e25 / w_bulk[2],
                rates_count.getRate(mc::ATTACHMENT) * 2.4e25 / w_bulk[2]
            );

        } else {
            std::printf(
                "\n Collisions: %lu\n"
                " Electrons: %i\n"
                " Cations: %i\n"
                " Anions: %i\n"
                " Mean Energy: %.2e eV\n"
                " Current Time: %.3e ns\n",

                collisions,
                mean.back().get_particles()[mc::ELECTRONS],
                mean.back().get_particles()[mc::CATIONS],
                mean.back().get_particles()[mc::ANIONS],
                mean.back().get_energy(),
                t.back() * 1e9

            );
        }
    }
}

void MonteCarlo::saveResults(const int64_t duration) const {
    // Prints the final results of the simulation to a file

    // Generate file .txt to save results
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::tm* local_time = std::localtime(&time_t);
    std::stringstream ss;
    ss << "results/" << std::put_time(local_time, "%Y-%m-%d_%H-%M") << ".txt";
    std::string filename = ss.str();
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not create results file: " << filename << std::endl;
        return;
    }
    
    file << std::fixed << std::setprecision(6);
    
    // Header
    file << "# Monte Carlo Electron Transport Simulation Results\n";
    file << "# Generated on: " << std::ctime(&time_t);
    file << "# ================================================\n\n";
    
    // Simulation parameters
    file << "[SIMULATION_PARAMETERS]\n";
    file << "EN = " << EN << " Td\n";
    file << "N0 = " << N0 << "\n";
    file << "N = " << N << " m^-3\n";
    file << "col_max = " << col_max << "\n";
    file << "col_equ = " << col_equ << "\n";
    file << "conserve = " << (conserve ? "true" : "false") << "\n";
    file << "isotropic = " << (isotropic ? "true" : "false") << "\n";
    file << "W = " << W << "\n";
    file << "E_max = " << E_max << " eV\n\n";
    
    // Gas composition
    file << "[GAS_COMPOSITION]\n";
    for (size_t i = 0; i < gas.size(); i++) {
        file << "gas[" << i << "] = " << gas[i] << "\n";
        file << "mix[" << i << "] = " << mix[i] << "\n";
    }
    file << "\n";
    
    // Energy data
    file << "[ENERGY_DATA]\n";
    file << "E_mean = " << E.get_E_mean() << " eV\n\n";
    
    // Bulk transport data
    if (!bulk.is_empty()) {
        const auto& w_bulk = bulk.get_w();
        const auto& w_bulk_err = bulk.get_w_err();
        const auto& DN_bulk = bulk.get_DN();
        const auto& DN_bulk_err = bulk.get_DN_err();
        
        file << "[BULK_TRANSPORT_DATA]\n";
        file << "w_x = " << w_bulk[0] << " m/s\n";
        file << "w_y = " << w_bulk[1] << " m/s\n";
        file << "w_z = " << w_bulk[2] << " m/s\n";
        file << "w_x_err = " << w_bulk_err[0] << " m/s\n";
        file << "w_y_err = " << w_bulk_err[1] << " m/s\n";
        file << "w_z_err = " << w_bulk_err[2] << " m/s\n";
        file << "DN_x = " << DN_bulk[0] << " m^2/s\n";
        file << "DN_y = " << DN_bulk[1] << " m^2/s\n";
        file << "DN_z = " << DN_bulk[2] << " m^2/s\n";
        file << "DN_x_err = " << DN_bulk_err[0] << " m^2/s\n";
        file << "DN_y_err = " << DN_bulk_err[1] << " m^2/s\n";
        file << "DN_z_err = " << DN_bulk_err[2] << " m^2/s\n\n";
    }
    
    // Flux transport data
    const auto& w_flux = flux.get_w();
    const auto& DN_flux = flux.get_DN();
    
    file << "[FLUX_TRANSPORT_DATA]\n";
    file << "w_x = " << w_flux[0] << " m/s\n";
    file << "w_y = " << w_flux[1] << " m/s\n";
    file << "w_z = " << w_flux[2] << " m/s\n";
    file << "DN_x = " << DN_flux[0] << " m^2/s\n";
    file << "DN_y = " << DN_flux[1] << " m^2/s\n";
    file << "DN_z = " << DN_flux[2] << " m^2/s\n\n";

    // Transport coefficients
    if (!bulk.is_empty()) {
        const auto& w_bulk = bulk.get_w();
        file << "[TRANSPORT_COEFFICIENTS]\n";
        file << "alpha = " << rates_count.getRate(mc::IONIZATION) * 2.4e25 / w_bulk[2] << " m^-1\n";
        file << "eta = " << rates_count.getRate(mc::ATTACHMENT) * 2.4e25 / w_bulk[2] << " m^-1\n\n";
    }
    
    // Reaction rates
    file << "[REACTION_RATES]\n";
    file << "# Counted rates:\n";
    file << "effective_count = " << rates_count.getRate(mc::EFFECTIVE) * 1e12 << "e+12 m^3/s\n";
    file << "ionization_count = " << rates_count.getRate(mc::IONIZATION) * 1e12 << "e+12 m^3/s\n";
    file << "attachment_count = " << rates_count.getRate(mc::ATTACHMENT) * 1e12 << "e+12 m^3/s\n";
    file << "# Rates computed by convolution:\n";
    file << "effective_conv = " << rates_conv.getRate(mc::EFFECTIVE) * 1e12 << "e+12 m^3/s\n";
    file << "ionization_conv = " << rates_conv.getRate(mc::IONIZATION) * 1e12 << "e+12 m^3/s\n";
    file << "attachment_conv = " << rates_conv.getRate(mc::ATTACHMENT) * 1e12 << "e+12 m^3/s\n\n";

    // Final state
    file << "[FINAL_STATE]\n";
    file << "iterations = " << t.size() << "\n";
    file << "collisions = " << collisions << "\n";
    file << "electrons = " << v.size() << "\n";
    file << "cations = " << r[mc::CATIONS].size() << "\n";
    file << "anions = " << r[mc::ANIONS].size() << "\n";
    file << "steady_state_time = " << T_sst * 1e9 << " ns\n";
    file << "final_time = " << t.back() * 1e9 << " ns\n";
    file << "convergence_status = " << converge << "\n";
    if (duration >= 3600) {
        const int hours = duration / 3600;
        const int minutes = (duration % 3600) / 60;
        file << "convergence_time = " << hours << " hours " << minutes << " minutes\n\n";
    } else {
        const int minutes = duration / 60;
        file << "convergence_time = " << minutes << " minutes\n\n";
    }
    
    // Energy losses
    file << "[ENERGY_LOSSES]\n";
    file << "elastic = " << EnergyLossElastic << " eV\n";
    file << "inelastic = " << EnergyLossInelastic << " eV\n";
    file << "ionization = " << EnergyLossIonization << " eV\n\n";

    // Energy distribution
    const auto& energy_grid = E.get_energy();
    const auto& EEDF = E.get_EEDF();
    const auto& EEPF = E.get_EEPF();
    file << "[ENERGY_DISTRIBUTION]\n";
    file << "# energy[eV] EEDF EEPF\n";
    //size_t n_points = std::min(static_cast<size_t>(100), energy_grid.size());
    for (size_t i = 0; i < energy_grid.size(); i++) {
        file << energy_grid[i] << " " << EEDF[i] << " " << EEPF[i] << "\n";
    }
    file << "\n";
    
    file.close();
    std::cout << " Results saved to: " << filename << "\n" << std::endl;
}