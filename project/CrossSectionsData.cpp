#include "CrossSectionsData.hpp"

CrossSectionsData::CrossSectionsData(const std::vector<std::string> & species, const double E_max, 
                                     const std::vector<double> & mix, const double N):
    gas(species), n_react(species.size(),0), nu_max(0.0)
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

    // Make sure that the energy vector contains the endpoints of the energy grid of the simulation:
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
        import_Xsec_data(offset, gas[i]);
        offset += n_react[i];   // Update the offset for the next specie
    }

    // Compute the maximal collision frequency:
    maximalCollFreq(mix, N);
}

// Counts the number of reactions for each specie and the number of cross-sections data points
// and fills the energy vector with the energy levels from each reaction, sorted and deduplicated.
void CrossSectionsData::count_react_and_fill_energy(const size_t i){
    // Path to the cross-section data file
    std::string path = "./Xsec/" + gas[i] + "/" + gas[i] + ".txt";

    std::ifstream file;               
    std::string line;     
    bool flag=false;      
    unsigned counter=0;       

    file.open(path,std::ifstream::in);
    if(file.is_open()){
        while(getline(file,line))
        {
            if (line.compare("EFFECTIVE")==0 || line.compare("IONIZATION")==0 ||
                line.compare("ATTACHMENT")==0 || line.compare("EXCITATION")==0 ||
                line.compare("ELASTIC") == 0){
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
                    std::stable_sort(energy.begin(), energy.end());
                    energy.erase(std::unique(energy.begin(), energy.end(),
                        [](double a, double b){ return std::fabs(a - b) < 1e-20; }),
                        energy.end());
                    // Note: energy values are in the range (0 eV - 40 eV) approximately,
                    // so energy levels that differ by less than 0.1 eV are considered equal.
                }
            }
        }
    }
}

// void CrossSectionsData::import_Xsec_data(const size_t i){
void CrossSectionsData::import_Xsec_data(const size_t offset, const std::string & current_specie){

    // Path to the cross-section data file
    std::string path = "./Xsec/" + current_specie + "/" + current_specie + ".txt";

    // Map interactions indexes with their actual name in .txt data files
    const std::map<std::string, mc::InteractionType> int_map = {
        {"EFFECTIVE", mc::EFFECTIVE},
        {"IONIZATION", mc::IONIZATION},
        {"ATTACHMENT", mc::ATTACHMENT},
        {"EXCITATION", mc::EXCITATION},
        {"ELASTIC", mc::ELASTIC}
    };

    std::ifstream file;               
    std::string line;                  
    std::string avg_energy,formula;
    bool flag = true;
    unsigned counter=0;
    std::vector<double> xx, yy;      // energy levels and cross-section values

    // Import data relative to the current specie in Xsec_tool[i]: 
    file.open(path,std::ifstream::in);
    if(file.is_open()){

        int j=0;                         // reaction index

        while(getline(file,line)) {
            // Set interaction type and average energy:
            if (line.compare("EFFECTIVE")==0 || line.compare("IONIZATION")==0 ||
                line.compare("ATTACHMENT")==0 || line.compare("EXCITATION")==0 ||
                line.compare("ELASTIC") == 0){

                flag = true;
                auto it = int_map.find(line);

                // Set interaction type
                //Xsections[i][j].interact = it->second;
                Xsections[offset + j].interact = it->second;
                
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
    }
}

// Linear interpolation function
void CrossSectionsData::linear_interpolation(std::vector<double>& x, std::vector<double>& y, std::vector<double>& result) {

    if (x.empty() || y.empty() || x.size() != y.size()) {
        throw std::invalid_argument("Input vectors x and y must be non-empty and of the same size.");
    }

   // Ensure that the energy vector contains the energy level E = 0.0:
   double E_max = energy.back();

    if( x[0] > 0.0) {
        x.insert(x.begin(), 0.0);
        y.insert(y.begin(), y.front());
    }
    if( x.back() < E_max) {
        x.push_back(E_max);
        y.push_back(y.back());
    }

    for (const double q : energy) {
        size_t i = 0;

        if (q <= x.front()) {
            i = 0;
        } else if (q >= x.back()) {
            i = x.size() - 2;
        } else {
            while (i < x.size() - 1 && q > x[i + 1]) i++;
        }

        double t = (q - x[i]) / (x[i + 1] - x[i]);
        result.push_back(y[i] + t * (y[i + 1] - y[i]));
    }
}

// Compute the maximal collision frequency:
void CrossSectionsData::maximalCollFreq(const std::vector<double> & mix, const double N) {

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