#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <regex>
#include <iomanip>
#include <map>
#include <limits>  // Per NAN

struct SimulationResult {
    double convergence_time_seconds;
    double E_mean;
    double w_z_bulk;
    double DN_z_bulk;
    double w_z_flux;
    double DN_z_flux;
    double ionization_count;
    double attachment_count;
    double ionization_conv;
    double attachment_conv;
    
    // Constructor con valori di default
    SimulationResult() : convergence_time_seconds(-1), E_mean(-1), w_z_bulk(-1), 
                        DN_z_bulk(-1), w_z_flux(-1), DN_z_flux(-1), 
                        ionization_count(-1), attachment_count(-1), 
                        ionization_conv(-1), attachment_conv(-1) {}
};

class MonteCarloAnalyzer {
private:
    struct Statistics {
        double mean;
        double variance;
        double std_dev;
        int count;
        double min_val;
        double max_val;
    };

public:
    // Converte il tempo da formato "X minutes, Y seconds" a secondi totali
    double parseConvergenceTime(const std::string& time_str) {
        std::regex time_regex(R"((\d+)\s+minutes?,\s+(\d+)\s+seconds?)");
        std::smatch matches;
        
        if (std::regex_search(time_str, matches, time_regex)) {
            int minutes = std::stoi(matches[1].str());
            int seconds = std::stoi(matches[2].str());
            return minutes * 60.0 + seconds;
        }
        
        // Prova solo secondi
        std::regex sec_regex(R"((\d+(?:\.\d+)?)\s+seconds?)");
        if (std::regex_search(time_str, matches, sec_regex)) {
            return std::stod(matches[1].str());
        }
        
        return -1; // Errore nel parsing
    }
    
    // Estrae un valore numerico da una riga del tipo "key = value unit"
    double extractValue(const std::string& line, const std::string& key) {
        size_t key_pos = line.find(key + " =");
        if (key_pos == std::string::npos) return -1;
        
        size_t equals_pos = line.find("=", key_pos);
        if (equals_pos == std::string::npos) return -1;
        
        std::string value_part = line.substr(equals_pos + 1);
        std::istringstream iss(value_part);
        double value;
        
        if (iss >> value) {
            return value;
        }
        
        return -1;
    }
    
    // Legge un file di risultati e estrae i valori di interesse
    SimulationResult readResultFile(const std::string& filename) {
        SimulationResult result;
        std::ifstream file(filename);
        std::string line;
        std::string current_section = "";
        
        if (!file.is_open()) {
            std::cerr << "Errore: impossibile aprire " << filename << std::endl;
            return result;
        }
        
        while (std::getline(file, line)) {
            // Rimuovi spazi iniziali e finali
            line.erase(0, line.find_first_not_of(" \t"));
            line.erase(line.find_last_not_of(" \t") + 1);
            
            // Identifica la sezione corrente
            if (line.find("[") == 0 && line.find("]") != std::string::npos) {
                current_section = line;
                continue;
            }
            
            // Salta commenti e righe vuote
            if (line.empty() || line[0] == '#') continue;
            
            // Parsing del convergence_time
            if (line.find("convergence_time =") != std::string::npos) {
                std::string time_part = line.substr(line.find("=") + 1);
                result.convergence_time_seconds = parseConvergenceTime(time_part);
            }
            
            // Energy data
            else if (line.find("E_mean =") != std::string::npos) {
                result.E_mean = extractValue(line, "E_mean");
            }
            
            // Bulk transport data
            else if (current_section == "[BULK_TRANSPORT_DATA]") {
                if (line.find("w_z =") != std::string::npos) {
                    result.w_z_bulk = extractValue(line, "w_z");
                }
                else if (line.find("DN_z =") != std::string::npos) {
                    result.DN_z_bulk = extractValue(line, "DN_z");
                }
            }
            
            // Flux transport data
            else if (current_section == "[FLUX_TRANSPORT_DATA]") {
                if (line.find("w_z =") != std::string::npos) {
                    result.w_z_flux = extractValue(line, "w_z");
                }
                else if (line.find("DN_z =") != std::string::npos) {
                    result.DN_z_flux = extractValue(line, "DN_z");
                }
            }
            
            // Reaction rates
            else if (current_section == "[REACTION_RATES]") {
                if (line.find("ionization_count =") != std::string::npos) {
                    result.ionization_count = extractValue(line, "ionization_count");
                }
                else if (line.find("attachment_count =") != std::string::npos) {
                    result.attachment_count = extractValue(line, "attachment_count");
                }
                else if (line.find("ionization_conv =") != std::string::npos) {
                    result.ionization_conv = extractValue(line, "ionization_conv");
                }
                else if (line.find("attachment_conv =") != std::string::npos) {
                    result.attachment_conv = extractValue(line, "attachment_conv");
                }
            }
        }
        
        return result;
    }
    
    // Calcola statistiche per un vettore di valori
    Statistics calculateStatistics(const std::vector<double>& values) {
        Statistics stats = {0.0, 0.0, 0.0, 0, 0.0, 0.0};
        
        if (values.empty()) return stats;
        
        stats.count = values.size();
        stats.min_val = values[0];
        stats.max_val = values[0];
        
        // Calcola media e trova min/max
        double sum = 0.0;
        for (double val : values) {
            sum += val;
            stats.min_val = std::min(stats.min_val, val);
            stats.max_val = std::max(stats.max_val, val);
        }
        stats.mean = sum / stats.count;
        
        // Calcola varianza
        double variance_sum = 0.0;
        for (double val : values) {
            variance_sum += (val - stats.mean) * (val - stats.mean);
        }
        stats.variance = variance_sum / stats.count;
        stats.std_dev = std::sqrt(stats.variance);
        
        return stats;
    }
    
    // Analizza tutti i file nella cartella results/
    void analyzeResults(const std::string& results_dir = "results/") {
        std::vector<SimulationResult> all_results;
        std::vector<std::string> filenames;
        
        // Leggi tutti i file result*.txt
        for (int i = 1; i <= 50; ++i) {
            std::string filename = results_dir + "result" + std::to_string(i) + ".txt";
            SimulationResult result = readResultFile(filename);
            
            // Verifica che almeno alcuni valori siano stati letti correttamente
            if (result.convergence_time_seconds > 0 || result.E_mean > 0) {
                all_results.push_back(result);
                filenames.push_back("result" + std::to_string(i) + ".txt");
                std::cout << "Letto " << filename << " - convergence_time: " 
                         << result.convergence_time_seconds << "s" << std::endl;
            } else {
                std::cout << "Attenzione: " << filename << " non contiene dati validi" << std::endl;
            }
        }
        
        if (all_results.empty()) {
            std::cerr << "Nessun file con dati validi trovato!" << std::endl;
            return;
        }
        
        std::cout << "\nTrovati " << all_results.size() << " file validi." << std::endl;
        
        // Estrai i valori per ogni variabile
        std::map<std::string, std::vector<double>> variables;
        
        for (const auto& result : all_results) {
            if (result.convergence_time_seconds > 0) 
                variables["convergence_time_s"].push_back(result.convergence_time_seconds);
            if (result.E_mean > 0) 
                variables["E_mean_eV"].push_back(result.E_mean);
            if (result.w_z_bulk > 0) 
                variables["w_z_bulk_ms"].push_back(result.w_z_bulk);
            if (result.DN_z_bulk > 0) 
                variables["DN_z_bulk_m2s"].push_back(result.DN_z_bulk);
            if (result.w_z_flux > 0) 
                variables["w_z_flux_ms"].push_back(result.w_z_flux);
            if (result.DN_z_flux > 0) 
                variables["DN_z_flux_m2s"].push_back(result.DN_z_flux);
            if (result.ionization_count > 0) 
                variables["ionization_count_m3s"].push_back(result.ionization_count);
            if (result.attachment_count > 0) 
                variables["attachment_count_m3s"].push_back(result.attachment_count);
            if (result.ionization_conv > 0) 
                variables["ionization_conv_m3s"].push_back(result.ionization_conv);
            if (result.attachment_conv > 0) 
                variables["attachment_conv_m3s"].push_back(result.attachment_conv);
        }
        
        // Crea il CSV con le statistiche
        std::ofstream csv_file("monte_carlo_statistics.csv");
        csv_file << std::scientific << std::setprecision(6);  // Notazione scientifica con 6 cifre significative
        
        csv_file << "Variable,Mean,Variance,StdDev,Min,Max,Count" << std::endl;
        
        for (const auto& [var_name, values] : variables) {
            if (!values.empty()) {
                Statistics stats = calculateStatistics(values);
                csv_file << var_name << ","
                        << stats.mean << ","
                        << stats.variance << ","
                        << stats.std_dev << ","
                        << stats.min_val << ","
                        << stats.max_val << ","
                        << stats.count << std::endl;
                
                std::cout << var_name << ": "
                         << "Media=" << stats.mean 
                         << ", StdDev=" << stats.std_dev
                         << ", Count=" << stats.count << std::endl;
            }
        }
        
        csv_file.close();
        std::cout << "\nStatistiche salvate in monte_carlo_statistics.csv" << std::endl;
        
        // Crea CSV con tutti i dati grezzi
        createRawDataCSV(all_results, filenames);
    }
    
private:
    void createRawDataCSV(const std::vector<SimulationResult>& all_results,
                         const std::vector<std::string>& filenames) {
        std::ofstream raw_csv("monte_carlo_raw_data.csv");
        raw_csv << std::scientific << std::setprecision(6);  // Notazione scientifica
        
        // Header
        raw_csv << "Simulation,convergence_time_s,E_mean_eV,w_z_bulk_ms,DN_z_bulk_m2s,"
                << "w_z_flux_ms,DN_z_flux_m2s,ionization_count_m3s,attachment_count_m3s,"
                << "ionization_conv_m3s,attachment_conv_m3s" << std::endl;
        
        // Dati
        for (size_t i = 0; i < all_results.size(); ++i) {
            const auto& result = all_results[i];
            raw_csv << filenames[i] << ",";
            
            // Uso operatore ternario con notazione scientifica diretta
            raw_csv << (result.convergence_time_seconds > 0 ? result.convergence_time_seconds : NAN) << ","
                   << (result.E_mean > 0 ? result.E_mean : NAN) << ","
                   << (result.w_z_bulk > 0 ? result.w_z_bulk : NAN) << ","
                   << (result.DN_z_bulk > 0 ? result.DN_z_bulk : NAN) << ","
                   << (result.w_z_flux > 0 ? result.w_z_flux : NAN) << ","
                   << (result.DN_z_flux > 0 ? result.DN_z_flux : NAN) << ","
                   << (result.ionization_count > 0 ? result.ionization_count : NAN) << ","
                   << (result.attachment_count > 0 ? result.attachment_count : NAN) << ","
                   << (result.ionization_conv > 0 ? result.ionization_conv : NAN) << ","
                   << (result.attachment_conv > 0 ? result.attachment_conv : NAN) << std::endl;
        }
        
        raw_csv.close();
        std::cout << "Dati grezzi salvati in monte_carlo_raw_data.csv" << std::endl;
    }
};

int main() {
    MonteCarloAnalyzer analyzer;
    
    std::cout << "Analisi dei risultati Monte Carlo..." << std::endl;
    std::cout << "Estraendo: convergence_time, E_mean, w_z (bulk/flux), DN_z (bulk/flux)," << std::endl;
    std::cout << "           ionization_count/conv, attachment_count/conv" << std::endl;
    std::cout << "Cercando file nella cartella results/" << std::endl;
    
    analyzer.analyzeResults();
    
    std::cout << "\nAnalisi completata!" << std::endl;
    std::cout << "File generati:" << std::endl;
    std::cout << "- monte_carlo_statistics.csv (statistiche per variabile)" << std::endl;
    std::cout << "- monte_carlo_raw_data.csv (tutti i dati per ulteriori analisi)" << std::endl;
    
    return 0;
}