//#include "utils/CrossSectionsData.hpp"
#include "../../include/utils/CrossSectionsData.hpp" // for testing

CrossSectionsData::CrossSectionsData(const std::vector<std::string> & species, const double E_max, 
                                     const std::vector<double> & mix, const double N):
    gas(species), n_react(species.size(),1), nu_max(0.0)
{
    // Check if the gas species are the same as the cross-section data
    std::array<std::string, 6> valid_species = {"H2", "H2O", "N", "N2", "O", "O2"};
    for(const auto & specie : gas) {
        if (std::find(std::begin(valid_species), std::end(valid_species), specie) == std::end(valid_species)) {
            throw std::invalid_argument("Invalid gas species: " + specie);
        }
    }

    const size_t n_species = gas.size();

    // Count the number of reactions for each specie and build a common energy grid:
    for( size_t i = 0; i < n_species; i++){
        count_react_and_fill_energy(i);
    }

    // Make sure that the energy range for XS data contains the energy levels of the simulation:
    if( energy[0] > 0.0) energy.insert(energy.begin(), 0.0);
    if( energy.back() < E_max) energy.push_back(E_max);
    const size_t en_sz = energy.size();

    // "offset" identifies the starting index of the current specie in Xsections:
    size_t offset = 0;
    // Fill every reaction with the cross section data:
    for( size_t i = 0; i < n_species; i++) {
        // Pre-allocate the memory for the Xsections vectors:

        for( size_t j = 0; j < n_react[i]; j++) {
            table tab;
            tab.specie_index = i; // Set specie index
            tab.section.reserve(en_sz);
            Xsections.push_back(tab);
        }

        // Fill Xsections by interpolating imported data over the energy grid:
        import_Xsec_data(offset, i);
        offset += n_react[i];   // Update the offset for the next specie
    }

    // Compute the maximal collision frequency:
    maximalCollFreq(mix, N);
}

const unsigned CrossSectionsData::get_n_react() const {
    unsigned sum = 0;
    for (auto n : n_react) sum += n;
    return sum - gas.size(); // Exclude effective xs data for each specie
}

void CrossSectionsData::count_react_and_fill_energy(const size_t i){
    // Counts the number of reactions for each specie and the number of cross-sections data points
    // and fills the energy vector with the energy levels from each reaction, sorted and deduplicated.

    // Path to the cross-section data file
    std::string path = "data/Xsec/" + gas[i] + "/" + gas[i] + ".txt";

    std::ifstream file;               
    std::string line;     
    bool flag=false;      
    unsigned counter=0;       

    file.open(path,std::ifstream::in);
    if(file.is_open()){
        while(getline(file,line))
        {
            if (line.compare("EFFECTIVE")==0 || line.compare("IONIZATION")==0 ||
                line.compare("ATTACHMENT")==0 || line.compare("EXCITATION")==0){
                n_react[i]++;
                flag = true;
            }
            
            if(flag){
                if(line[0]=='-') counter++;
                std::stringstream ss(line);
                double x, y;  // x = energy level, y = cross section value

                if(counter < 2) {
                    if(ss >> x >> y) {
                        energy.push_back(x);
                    }
                } else {
                    flag = false;
                    counter = 0;
                    // Sort and deduplicate the energy vector:
                    /*
                    std::stable_sort(energy.begin(), energy.end());
                    energy.erase(std::unique(energy.begin(), energy.end(),
                        [](double a, double b){ return std::fabs(a - b) < 1e-10; }),
                        energy.end());
                    */
                    std::set<double> s(energy.begin(),energy.end());    // sort + deduplicate using set
                    energy.assign(s.begin(),s.end());
                }
            }
        }
    }
}

void CrossSectionsData::import_Xsec_data(const size_t offset, const size_t specie_index){
    // Imports cross-section data from .txt files

    // Path to the cross-section data file
    std::string path = "data/Xsec/" + gas[specie_index] + "/" + gas[specie_index] + ".txt";

    // Map interactions indexes with their actual name in .txt data files
    const std::map<std::string, mc::InteractionType> int_map = {
        {"EFFECTIVE", mc::EFFECTIVE},
        {"IONIZATION", mc::IONIZATION},
        {"ATTACHMENT", mc::ATTACHMENT},
        {"EXCITATION", mc::EXCITATION}
    };

    std::ifstream file;               
    std::string line;                  
    std::string avg_energy,formula;
    bool flag = true;
    unsigned counter=0;
    std::vector<double> xx, yy;     // energy levels and cross-section values

    // Import data relative to the current specie in Xsec_tool[i]: 
    file.open(path,std::ifstream::in);
    if(file.is_open()){

        int j=0;                         // reaction index
        int j_effective = 0;             // index of effective cross-section

        while(getline(file,line)) {
            // Set interaction type and average energy:
            if (line.compare("EFFECTIVE")==0 || line.compare("IONIZATION")==0 ||
                line.compare("ATTACHMENT")==0 || line.compare("EXCITATION")==0){

                flag = true;
                auto it = int_map.find(line);
                
                // Set interaction type
                //Xsections[i][j].interact = it->second;
                Xsections[offset + j].interact = it->second;

                if( it->second == mc::EFFECTIVE) {
                    j_effective = j;
                }
                
                // Set reaction formula:
                getline(file,formula);
                //Xsections[i][j].react=formula;
                Xsections[offset + j].react=formula;
                
                // Set average energy:
                getline(file,avg_energy);
                std::stringstream ss_en(avg_energy);
                //ss_en >> Xsections[i][j].en_avg;
                ss_en >> Xsections[offset + j].en_avg;
            }


            // Set energy and cross sectional values:
            if(flag) {
                if(line[0]=='-') counter++;
                std::stringstream ss(line);
                double x,y;  // x = energy level; y = cross section value

                if(counter < 2) {
                    if(ss >> x >> y) {
                        // Store energy and cross-section values in the temporary vectors
                        xx.push_back(x);
                        yy.push_back(y);
                    }
                } else {
                    flag = false;
                    counter = 0;

                    // Interpolate cross-section values to the common energy grid
                    if( !xx.empty() && !yy.empty()) {
                        linear_interpolation(xx, yy, Xsections[offset + j].section);
                        // Clear temporary vectors for the next interaction
                        xx.clear();
                        yy.clear();
                        j++;
                    }
                }
            }
        }

        // Last, compute elastic cross sectional data (as in Vahedi et al.):
        Xsections[offset + j].interact = mc::ELASTIC;
        Xsections[offset + j].react = "Elastic";
        Xsections[offset + j].en_avg = 0.0;
        Xsections[offset + j].section.reserve(energy.size());

        const double beta = 1.0;  // (no correction)

        for(size_t i = 0; i < energy.size(); i++){
            // Compute energy-dependent correction coeff. (=1 for small energies)
            // const double E = energy[i];
            // beta = (E < 1e-6) ? 1.0 : 0.5 * ( E * std::log(1+E) ) / (E - std::log(1+E));
            Xsections[offset + j].section.push_back(beta * Xsections[offset + j_effective].section[i]);

            // Subtract excitation and ionization
            for(int k = 0; k < j; k++) {
                if(Xsections[offset + k].interact == mc::EXCITATION ||
                    Xsections[offset + k].interact == mc::IONIZATION)
                {
                    Xsections[offset + j].section[i] -= Xsections[offset + k].section[i];
                }
            }
            // avoid negative xs:
            Xsections[offset + j].section[i] = std::max(Xsections[offset+j].section[i],0.0);
        }
    }
}

void CrossSectionsData::linear_interpolation(std::vector<double>& x, std::vector<double>& y, std::vector<double>& result) {
    // Linear interpolation function

    // Ensure that the energy vector contains the energy level E = 0.0:
    // const double E_max = energy.back();
    const double x_min = x.front();
    const double x_max = x.back();
    const double y_min = y.front();
    const double y_max = y.back();

    /*
    // Make the extremes of "energy" coincide with the simulation's energy range
    if( x[0] > 0.0) {
        x.insert(x.begin(), 0.0);
        y.insert(y.begin(), y.front());
    }
    if( x.back() < E_max) {
        x.push_back(E_max);
        y.push_back(y.back());
    }
    */
    //size_t k = 0;
    for (const double& q : energy) {
        if (q <= x_min) {
            result.push_back(y_min);
        } else if (q >= x_max) {
            result.push_back(y_max);
        } else {
            auto it = std::upper_bound(x.cbegin(), x.cend(), q);
            const size_t k = std::min(
                static_cast<size_t>(it - x.cbegin() - 1),
                x.size() - 2
            );
            //while (k < x.size() - 1 && q > x[i + 1]) k++;
            const double t = (q - x[k]) / (x[k + 1] - x[k]);
            result.push_back(y[k] + t * (y[k + 1] - y[k]));
        }
    }
}

void CrossSectionsData::maximalCollFreq(const std::vector<double> & mix, const double N) {
    // Compute the maximal collision frequency:

    // Check if the mix vector is the same size as the gas vector
    if (mix.size() != gas.size()) {
        throw std::invalid_argument("Mix vector size must match the number of gas species");
    }

    // Iterate over all energy levels:
    for( size_t k = 0; k < energy.size(); k++) {
        double sigma_tot = 0.0;   // total cross-section for each energy level
        // Iterate over all reactions:
        for (table& tab : Xsections) {
            // Compute the total cross-section at k-esim energy level:
            if( tab.interact != mc::EFFECTIVE) {
                sigma_tot += mix[tab.specie_index] * tab.section[k];
            }
        }
        // Calculate the maximal collision frequency:
        double nu = N * sigma_tot * std::sqrt(2.0 * energy[k] * mc::q0 / mc::me); // s^-1
        if(nu > nu_max) nu_max = nu;
    }
    std::cout << "nu_max = " << nu_max << "\n" << std::endl; // -> nu_max è giusto
}

const std::vector<table> CrossSectionsData::get_Xsections( const size_t specie, const mc::InteractionType interaction) const {
    // Check if the specie and interaction indices are valid
    if (specie >= gas.size()) {
        throw std::out_of_range("Invalid specie index");
    }
    // Identify the offset for the current specie in Xsections:
    size_t offset = std::accumulate(n_react.begin(), n_react.begin() + specie, 0);

    // Return the cross-section data for the specified specie and interaction
    std::vector<table> result;
    for (size_t j = 0; j < n_react[specie]; j++) {
        if (Xsections[offset + j].interact == interaction) {
        result.push_back(Xsections[offset + j]);
        }
    }
    return result;
}

const std::vector<table> CrossSectionsData::get_Xsections( const std::string & specie, const std::string & interaction) const {
    size_t species_index = std::distance(gas.cbegin(), std::find(gas.cbegin(), gas.cend(), specie));

    // Check if the specie exists
    if (species_index == gas.size()) {
        throw std::invalid_argument("Invalid specie: " + specie);
    }

    // Map interaction strings to their corresponding enum values
    const std::map<std::string, mc::InteractionType> interaction_map = {
        {"EFFECTIVE", mc::EFFECTIVE},
        {"IONIZATION", mc::IONIZATION},
        {"ATTACHMENT", mc::ATTACHMENT},
        {"EXCITATION", mc::EXCITATION},
        {"ELASTIC", mc::ELASTIC}
    };

    // Check if the interaction exists
    auto it = interaction_map.find(interaction);
    if (it == interaction_map.end()) {
        throw std::invalid_argument("Invalid interaction: " + interaction);
    }
    // Get the cross-section data for the specified specie and interaction
    return get_Xsections(species_index, it->second);
}

void CrossSectionsData::remove_effective_xs() {
    // Within the MonteCarlo simulation, the full Xsections vector is accessed only by the classes
    // RateDataConv and CollisionData. In both cases, the presence effective cross section data is
    // not only unnecesary, but also it affects drastically the performance due to the need
    // of checking multiple times if the interaction is effective or not, before accessing the actual data.

    // Since the CrossSectionData class might have more general application as it only imports and stores
    // the cross-section data from .txt data files, this public method has been specifically defined
    // for the MonteCarlo collision-code application only and is called in the MonteCarlo class constructor. 

    // Remove effective xs data from the Xsections vector
    Xsections.erase(std::remove_if(Xsections.begin(), Xsections.end(),
        [](const table& t) { return t.interact == mc::EFFECTIVE; }), Xsections.end());
}