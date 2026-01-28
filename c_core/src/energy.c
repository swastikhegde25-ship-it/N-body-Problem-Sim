#include "energy.h"
#include "tree.h"
#include <math.h>
#include <omp.h>

// -------------------- ENERGY CALCULATIONS ------------------------
float calculate_potential_rec(Particle* p, int node_idx, SimConfig config) {
    Node* n = &node_pool[node_idx];
    float dx = n->x - p->x; 
    float dy = n->y - p->y; 
    float dz = n->z - p->z;
    float dist_sq = dx*dx + dy*dy + dz*dz;
    float dist = sqrtf(dist_sq + config.softening * config.softening);
    float width = n->max_x - n->min_x;
    float pot = 0.0f;

    if (width / dist < config.theta || n->particle_count == 1) {
        if (dist_sq > 0.0001f) {
             pot -= (config.G * n->mass) / dist;
        }
    } else {
        for(int i=0; i<8; i++) {
            if(n->children[i] != -1) {
                pot += calculate_potential_rec(p, n->children[i], config);
            }
        }
    }
    return pot;
}

EXPORT void get_energy_stats(Particle* particles, SimConfig config, EnergyStats* stats) {
    double total_ke = 0.0;
    double total_pe = 0.0;
    
    // Kinetic Energy
    int i;
    #pragma omp parallel for reduction(+:total_ke)
    for (i = 0; i < config.particle_count; i++) {
        float v2 = particles[i].vx * particles[i].vx + 
                   particles[i].vy * particles[i].vy + 
                   particles[i].vz * particles[i].vz;
        total_ke += 0.5f * particles[i].mass * v2;
    }

    // Potential Energy (Using Barnes-Hut Tree)
    if (global_root != -1) {
        #pragma omp parallel for reduction(+:total_pe)
        for (i = 0; i < config.particle_count; i++) {
            float phi = calculate_potential_rec(&particles[i], global_root, config);
            total_pe += 0.5f * particles[i].mass * phi; // 0.5 to avoid double counting
        }
    }

    stats->kinetic = total_ke;
    stats->potential = total_pe;
    stats->total_energy = total_ke + total_pe;
}