#pragma once
#include "common.h"

float calculate_potential_rec(Particle* p, int node_idx, SimConfig config);
EXPORT void get_energy_stats(Particle* particles, SimConfig config, EnergyStats* stats);