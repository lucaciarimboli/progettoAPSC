#ifndef CONFIG_PARSER_HPP
#define CONFIG_PARSER_HPP

#include <nlohmann/json.hpp>
#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <iostream>
#include <cmath>

class ConfigParser {
public:
    // Simulation parameters
    std::vector<std::string> gas_species;
    std::vector<double> gas_mixture;
    double EN_field;
    double pressure;
    double temperature;
    unsigned long initial_electrons;
    double energy_sharing;
    std::array<double, 3> initial_position;
    std::array<double, 3> initial_broadening;

    // Grid parameters
    double max_energy;
    double energy_step;

    // Tolerances and limits
    double drift_velocity_error;
    double diffusion_error;
    unsigned long min_collisions;
    unsigned long max_collisions;
    unsigned long max_electrons;

    // Flags
    bool conserve_electrons;
    bool isotropic_scattering;
    bool save_results;

    // Methods
    bool loadFromFile(const std::string& filename) {
        try {
            std::ifstream file(filename);
            if (!file.is_open()) {
                std::cerr << "Cannot open config file: " << filename << std::endl;
                return false;
            }

            nlohmann::json j;
            file >> j;

            // Load simulation parameters
            if (j.contains("simulation")) {
                const auto& sim = j["simulation"];
                gas_species = sim.value("gas_species", gas_species);
                gas_mixture = sim.value("gas_mixture", gas_mixture);
                EN_field = sim.value("EN_field", EN_field);
                pressure = sim.value("pressure", pressure);
                temperature = sim.value("temperature", temperature);
                initial_electrons = sim.value("initial_electrons", initial_electrons);
                energy_sharing = sim.value("energy_sharing", energy_sharing);
                
                // Load arrays
                if (sim.contains("initial_position")) {
                    auto pos = sim["initial_position"].get<std::vector<double>>();
                    if (pos.size() == 3) {
                        std::copy(pos.begin(), pos.end(), initial_position.begin());
                    }
                }
                
                if (sim.contains("initial_broadening")) {
                    auto std_pos = sim["initial_broadening"].get<std::vector<double>>();
                    if (std_pos.size() == 3) {
                        std::copy(std_pos.begin(), std_pos.end(), initial_broadening.begin());
                    }
                }
            }

            // Load grid parameters
            if (j.contains("grid")) {
                const auto& grid = j["grid"];
                max_energy = grid.value("max_energy", max_energy);
                energy_step = grid.value("energy_step", energy_step);
            }

            // Load tolerances and limits
            if (j.contains("tolerances")) {
                const auto& tol = j["tolerances"];
                drift_velocity_error = tol.value("drift_velocity_error", drift_velocity_error);
                diffusion_error = tol.value("diffusion_error", diffusion_error);
                min_collisions = tol.value("min_collisions", min_collisions);
                max_collisions = tol.value("max_collisions", max_collisions);
                max_electrons = tol.value("max_electrons", max_electrons);
            }

            // Load flags
            if (j.contains("flags")) {
                const auto& flags = j["flags"];
                conserve_electrons = flags.value("conserve_electrons", conserve_electrons);
                isotropic_scattering = flags.value("isotropic_scattering", isotropic_scattering);
                save_results = flags.value("save_results", save_results);
            }

            return true;

        } catch (const std::exception& e) {
            std::cerr << "Error parsing JSON: " << e.what() << std::endl;
            return false;
        }
    }
};

#endif // CONFIG_PARSER_HPP