#include "physics_api.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* used github copilot inline code suggestion to speed up development 
but all the code logic, structure, and implementation is my own work.
Took help for some parts of the logic from stackoverflow */

// --- Main Loop ---
EXPORT void init_simulation() {}

EXPORT void step_simulation(Particle* particles, SimConfig config) {
    int i; 

    // Build Tree (Global Sync at t=0) ---
    reset_tree();
    #pragma omp parallel for private(i)
    for (i = 0; i < config.particle_count; i++) {
        particles[i].morton_index = morton_3d(particles[i].x, particles[i].y, particles[i].z, -config.world_size, config.world_size);
    }
    qsort(particles, config.particle_count, sizeof(Particle), compare_particles);

    // Update the global root variable so energy function can access it
    global_root = alloc_node();
    node_pool[global_root].min_x = -config.world_size; node_pool[global_root].max_x = config.world_size;
    node_pool[global_root].min_y = -config.world_size; node_pool[global_root].max_y = config.world_size;
    node_pool[global_root].min_z = -config.world_size; node_pool[global_root].max_z = config.world_size;

    for (i = 0; i < config.particle_count; i++) insert_node(global_root, i, &particles[i]);
    finalize_nodes(global_root);

    // Assign Timesteps ---
    assign_timesteps(particles, config);

    // Block Substep Execution ---
    int max_substeps = 1 << config.max_level;

    for (int s = 0; s < max_substeps; s++) {
        #pragma omp parallel for private(i)
        for (i = 0; i < config.particle_count; i++) {
            // Level L runs every (max_substeps / 2^L) ticks
            int interval = 1 << (config.max_level - particles[i].level);
            
            if (s % interval == 0) {
                leapfrog_active_step(particles, i, particles[i].current_dt, global_root, config);
            }
        }
    }
}