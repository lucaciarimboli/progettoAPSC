#include "cross_s.hpp"
#include "MonteCarlo.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <cmath>
#include <utility>
#include <random>

enum ParticleType {ELECTRONS = 0, CATIONS, ANIONS, PARTICLES_TYPES};
typedef std::vector<std::array<double,3>> POSITIONS;

std::array<POSITIONS,PARTICLES_TYPES> r;

// array of initial mean position of initial gaussian distributed electrons in x,y and z direction
std::array<double,3> pos_xyz  = {0, 0, 0};
// array of initial broadening of initial gaussian distributed electrons in x,y and z direction
std::array<double,3> sigma_xyz = {0, 0, 0};


