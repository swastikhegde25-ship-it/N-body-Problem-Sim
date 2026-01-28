#pragma once
#include "common.h"

void calculate_force(Particle* p, int node_idx, SimConfig config, float* fx, float* fy, float* fz);
void resolve_collisions(Particle* particles, int p_idx, int node_idx);
void assign_timesteps(Particle* particles, SimConfig config);
void leapfrog_active_step(Particle* particles, int p_idx, float dt, int root_idx, SimConfig config);