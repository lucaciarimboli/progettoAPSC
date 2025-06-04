#include "MonteCarlo.hpp"

void MonteCarlo::checkFractionSum(){

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
    MolMass mol_mass;
    for(const std::string & sub : gas){
        mol_mass.set_substance(sub);
        mol_mass.Compute_M();
        const auto& M_vec = mol_mass.get_M();
        for (size_t i = 0; i < M_vec.size(); ++i) {
            mgas[i] = M_vec[i] / (mc::Na * 1000);
        }
    }
}

std::pair<double,double> MonteCarlo::velocity2energy(const std::array<double,3> & v){
    // calculates absolute value of velocity abs_v and energy
    // E_in_eV in eV for one electron.
    double abs_v = std::sqrt(std::inner_product(v.cbegin(), v.cend(), v.cbegin(), 0.0));
    double E_in_eV =  0.5 * mc::me *  abs_v * abs_v / mc::q0;
    return std::make_pair(abs_v,E_in_eV);
}

void MonteCarlo::initialParticles(const std::array<double,3> & pos_xyz, const std::array<double,3> & sigma_xyz){

    // Initialize mean values
    MeanData m(pos_xyz,sigma_xyz,N0);
    
    mean.clear();
    mean.emplace_back(m);

    // Allocate memory to electron population positions:
    r[mc::ELECTRONS].resize(N0);

    // Initialize cations and anions populations as empty:
    r[mc::CATIONS].clear();
    r[mc::ANIONS].clear();

    // Random number generators
    std::default_random_engine gen;
    std::normal_distribution<double> randn(0.0,1.0);    // standard Gaussian distribution
    
    // Initialize electrons position
    for( auto &re : r[mc::ELECTRONS]){
        for( size_t i=0; i<=2; i++){ 
            // Initialize to pos_xyz + noise if sigma_xyz != 0
            re[i] = pos_xyz[i] + sigma_xyz[i] * randn(gen); 
        } 
    }

    // Electrons velocities are initialized with 0.0, add noise:
    // std::uniform_real_distribution<> randu(0.0,1.0); // [0,1] uniform distribution
    // for( auto &vel : v){
    //     for( size_t i=0; i<=2; i++){            
             // vel[i] += 1e-6 * randu(gen);
    //     } 
    // }
}

void MonteCarlo::set_E(const double EN){
    // Works for the case: EN constant, uniform and user defined.
    // "solvePoisson()" is an extension of this method.
    E_x = 0;
    E_y = 0;
    E_z = EN * N * 1e-21; // set E // z direction !
}

double MonteCarlo::random(){
        // Generates a random number from a U[0,1]
        std::default_random_engine gen;
        std::uniform_real_distribution<> randu(0.0,1.0); // [0,1] uniform distribution
        return randu(gen);
}

void MonteCarlo::freeFlight(){
    // performs non-collissional flight for electrons in electric field
    dt = - std::log(random()) / Xsec.get_nu_max(); // generates time step

    // Update vector time:
    t.emplace_back(t.back() + dt);
    
    double ne = v.size(); // number of electrons
    double dt2 = std::pow(dt,2); // to avoid redundant calculations 

    if(T_sst > 0.0){
        count_sst++;
        v_int.resize(ne);
        v2_int.resize(ne);
        double dt3 = std::pow(dt,3);
        for(size_t i = 0; i < ne; i++){     // i indicates the electron
            for(size_t j = 0; j < 3; j++){  // j indicates the coordinate (x,y,z)
                // integrated velocity
                v_int[i][j] = v[i][j]*dt + 0.5 * a[j] * dt2;
                // integrated velocity squared
                v2_int[i][j] = std::pow(v[i][j],2)*dt + a[j]*v[i][j] * dt2 + std::pow(a[j],2)/3 * dt3;
            }
        }
    }

    // Update space and velocity:
    for(size_t i = 0; i < ne; i++){
            for(size_t j = 0; j < 3; j++){
                r[mc::ELECTRONS][i][j] += v[i][j]*dt + 0.5 * a[j] * dt2;
                v[i][j] += a[j] * dt;
        }
    }

    //counter++;
}

void MonteCarlo::collectMeanData(){
    // Update mean vector for the new time step
    mean.emplace_back(mean.back().get_particles(), r[mc::ELECTRONS], v);
}

void MonteCarlo::updateEnergyData(){
    size_t ne = v.size();

    t_total += dt * ne; // sum of all times for all electrons
    std::vector<double> E_in_eV;
    E_in_eV.reserve(ne);
    std::transform(v.begin(), v.end(), std::back_inserter(E_in_eV), [this](const std::array<double, 3>& vi) {
        return velocity2energy(vi).second;
    });

    E.update_energy(E_in_eV,dt,ne,t_total);
    E.compute_distribution_function();  
}

void MonteCarlo::updateFluxData(){
    flux.compute_drift_velocity(v_int,t_total);
    flux.compute_diffusion_const(r[mc::ELECTRONS],v,N);
}

void MonteCarlo::updateBulkData(){
    bulk.update_bulk(t,count_sst,mean,N);
}

void MonteCarlo::updateReactionRates(){

    // 1. REACTION RATES BY COUNTING:
    rates_count.setTime(t,count_sst);
    rates_count.setParticles(mean,count_sst);
    rates_count.computeRates();

    // 2. REACTION RATES BY CONVOLUTION:
    rates_conv.setEnergy(E);
    rates_conv.computeRates();
};

void MonteCarlo::updateCollisionMatrix(){
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
    
    // Build collision matrix and compute indeces:
    //C.ComputeIndeces(num_particles, Xsec, E_in_eV, v_abs, mix, N);
    C.ComputeIndeces(num_particles, Xsec, E_in_eV, mix, N); // I express v_abs in terms of E_in_eV inside the computations to avoid allocating memory for it
    // Update total number of collisions:
    collisions += C.getCollisions();
}

void MonteCarlo::performCollision(const std::string & type){
    // Perform collisions according to the type of collision
    const std::vector<size_t> & ind = C.get_ind(type);
    if(ind.empty()) return; // in case no collisions occoured this time step

    const std::vector<double> & Mass = C.getMass();
    const std::vector<double> & Loss = C.getLoss();

    if(type == "ELASTIC"){
        elasticCollision(ind, Mass);
    }
    else if(type == "EXCITATION"){
        inelasticCollision(ind, Loss);
    }
    else if(type == "IONIZATION"){
        ionizationCollision(ind, Loss);
    }
    else if(type == "ATTACHMENT"){
        attachmentCollision(ind);
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
        double phi = 2 * M_PI * random();
        sin_phi = std::sin(phi);
        cos_phi = std::cos(phi);

        // Randomly generate xsi: electron scattering angle 
        if(isotropic) cos_xsi = 1 - 2 * random();
        else cos_xsi = (2 + v2e.second - 2 * (1+v2e.second) * random()) / v2e.second;
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
        E_2 += std::max(0.0, v2e.second*(1 - 2*mc::me/Mass[el_index] * (1 - cos_xsi)));

        // Update velocity:
        double v2_abs = std::sqrt(2.0 * E_2 * mc::q0 / mc::me);
        for (int j = 0; j < 3; j++) {
            v[el_index][j] = v2_abs * e_2[j];
        }
    }

    // Update total energy loss counter:
    EnergyLossElastic += E_1 - E_2;
}

void MonteCarlo::inelasticCollision(const std::vector<size_t> & ind, const std::vector<double> & Loss){

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
        double phi = 2 * M_PI * random();
        sin_phi = std::sin(phi);
        cos_phi = std::cos(phi);

        // Randomly generate xsi: electron scattering angle 
        if(isotropic) cos_xsi = 1 - 2 * random();
        else cos_xsi = (2 + v2e.second - 2 * (1+v2e.second) * random()) / v2e.second;
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
        E_2 += std::max(0.0, v2e.second - Loss[el_index]);

        // Update velocity:
        double v2_abs = std::sqrt(2.0 * E_2 * mc::q0 / mc::me);
        for (int j = 0; j < 3; j++) {
            v[el_index][j] = v2_abs * e_2[j];
        }
    }    

    // Update total energy loss counter:
    EnergyLossInelastic += E_1 - E_2;
}

void MonteCarlo::ionizationCollision(const std::vector<size_t> & ind, const std::vector<double> & Loss){

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
        double phi = 2 * M_PI * random();
        sin_phi = std::sin(phi);
        cos_phi = std::cos(phi);

        // Randomly generate xsi: electron scattering angle 
        if(isotropic) cos_xsi = 1 - 2 * random();
        else cos_xsi = (2 + v2e.second - 2 * (1+v2e.second) * random()) / v2e.second;
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
        E_2 += std::max(0.0, v2e.second - Loss[el_index]);

        // Update velocity:
        double v2_abs = std::sqrt(2.0 * W * E_2 * mc::q0 / mc::me); // W is the energy sharing in ionizing collision
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
            int random_index = static_cast<int>(random() * (r[mc::ELECTRONS].size()));
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
            int random_index = static_cast<int>(random() * (r[mc::ELECTRONS].size()));
            // clone it to compensate the removed one:
            r[mc::ELECTRONS].push_back(r[mc::ELECTRONS][random_index]);
            v.push_back(v[random_index]);
        }
    }

    // Update number of particles in mean data:
    mean.back().add_new_particles({- delta_Ne, 0, delta_Ne});
}

void MonteCarlo::checkSteadyState(){

    // If sst has not been reached yet, check after last iteration:
    if(count_sst == 0 && collisions/1e6 >= line && collisions >= col_equ){

        // check if the interval 80-90 % of energy data is larger than the interval 90-100%:
        int N = mean.size();
        size_t n = std::round( N / 10.0);

        double sum1 = 0.0, sum2 = 0.0;
        for (size_t i = N - 2 * n; i < N - n; i++) {
            sum1 += mean[i].get_energy();
        }
        for (size_t i = N - n; i < N; i++) {
            sum2 += mean[i].get_energy();
        }

        if(sum1 >= sum2){
            // Steady state has been reached
            T_sst = t.back();
            std::cout << "Steady state reached at t = " << T_sst << " ms\n";

            // For debugging purposes:
            std::cout << "Number of iterations until steady state: " << iii << "\n";
            std::cout << "Number of collisions until steady state: " << collisions << "\n";
            std::cout << "Number of electrons at steady state: " << v.size() << "\n";

            //counter = 0;
            collisions = 0;
            line = 1;
        }
    }
}

bool MonteCarlo::endSimulation() {
    // End simulation if too many electrons
    if (!conserve && v.size() > Ne_max) {
        converge = 1;
        std::cout << "Simulation ended: maximum number of electrons reached\n";
        return true;
    }

    // End simulation if no electrons
    if (!conserve && v.empty()) {
        converge = 2;
        std::cout << "Simulation ended: number of electrons is zero\n";
        return true;
    }

    // End simulation if relative errors in w and D are below thresholds
    if (!(bulk.is_empty())) {
        const std::array<double, 3> & w_bulk = bulk.get_w();
        const std::array<double, 3> & w_bulk_err = bulk.get_w_err();
        const std::array<double, 3> & DN_bulk = bulk.get_DN();
        const std::array<double, 3> & DN_bulk_err = bulk.get_DN_err();  
        if (std::abs(w_bulk_err[2] / w_bulk[2]) < w_err &&
            std::abs(DN_bulk_err[2] / DN_bulk[2]) < DN_err) {

            converge = 0;
            std::cout << "Simulation ended: errors in w < " 
                      << 100 * w_err << "% and D < "
                      << 100 * DN_err << "%\n";
            return true;
        }
    }

    // End simulation if number of collisions exceeds maximum
    if (collisions > col_max) {
        converge = 3;
        std::cout << "Simulation ended: maximum number of collisions reached\n";
        return true;
    }

    return false;

    // For debugging purposes:
    iii++;
}

void MonteCarlo::printOnScreen() {
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
                " Werr: %i"
                " DNerr %i"
                " collisions: %i"
                " electrons: %i"
                " E: %.3e eV"
                " w_bulk: %.3e m/s"
                " w_flux: %.3e m/s"
                " DN_bulk: %.2e (ms)^-1"
                " DN_flux: %.2e (ms)^-1"
                " Reff_count: %.2e m^3/s"
                " Reff_calc: %.2e m^3/s"
                " Alpha: %.3e m^-1"
                " Eta: %.3e m^-1\n",

                static_cast<int>(std::abs(w_bulk_err[2] / w_bulk[2])),
                static_cast<int>(std::abs(DN_bulk_err[2] / DN_bulk[2])),
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
                " collisions: %i"
                " electrons: %i"
                " mean energy: %.2e\n"


                // for debugging purposes:
                " ITERATION NUMBER: %i\n"
                " electrons according to vel: %i"
                " electrons according to pos: %i\n"
                " number of cations (pos): %i"
                " number of cations (mean): %i\n"
                " number of anions (pos): %i "
                " number of anions (mean): %i\n"
                " current time: %.3e ms"
                " dt: %.3e ms\n"
                " count_sst: %i \n\n",

                collisions,
                mean.back().get_particles()[mc::ELECTRONS],
                mean.back().get_energy(),

                // For debugging purposes:
                iii,
                v.size(),
                r[mc::ELECTRONS].size(),
                r[mc::CATIONS].size(),
                mean.back().get_particles()[mc::CATIONS],
                r[mc::ANIONS].size(),
                mean.back().get_particles()[mc::ANIONS],
                t.back(),
                dt,
                count_sst

            );
        }
    }
}
