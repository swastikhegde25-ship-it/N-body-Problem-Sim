#include "engine.h"
#include "tree.h"
#include <math.h>
#include <omp.h>

void calculate_force(Particle* p, int node_idx, SimConfig config, float* fx, float* fy, float* fz) {
    Node* n = &node_pool[node_idx];
    float dx = n->x - p->x; float dy = n->y - p->y; float dz = n->z - p->z;
    float dist_sq = dx*dx + dy*dy + dz*dz;
    float dist = sqrtf(dist_sq + config.softening * config.softening);
    float width = n->max_x - n->min_x;

    // Condition to use this node as an approximation
    if (width / dist < config.theta || n->particle_count == 1) {
        if (dist_sq > 0.0001f) {
             float f = (config.G * p->mass * n->mass) / (dist_sq * dist);
             *fx += f * dx; *fy += f * dy; *fz += f * dz;
        }
    } else {
        // Must recurse. 
        int recursed_any = 0;
        for(int i=0; i<8; i++) {
            if(n->children[i] != -1) {
                calculate_force(p, n->children[i], config, fx, fy, fz);
                recursed_any = 1;
            }
        }
         
        // accumulated mass. Otherwise, this mass is ignored (gravity disappears).
        if (!recursed_any && dist_sq > 0.0001f) {
             float f = (config.G * p->mass * n->mass) / (dist_sq * dist);
             *fx += f * dx; *fy += f * dy; *fz += f * dz;
        }
    }
}

void resolve_collisions(Particle* particles, int p_idx, int node_idx) {
    Node* n = &node_pool[node_idx];
    Particle* p = &particles[p_idx];

    if (n->particle_count == 1) {
        int other = n->first_particle;
        if (other != -1 && other != p_idx) {
            float pdx = particles[other].x - p->x;
            float pdy = particles[other].y - p->y;
            float pdz = particles[other].z - p->z;
            float dist_sq = pdx*pdx + pdy*pdy + pdz*pdz;
            
            if (dist_sq < 16.0f && dist_sq > 0.0001f) {
                float dist = sqrtf(dist_sq);
                float overlap = 0.5f * (4.0f - dist);
                float nx = pdx / dist; float ny = pdy / dist; float nz = pdz / dist;
                p->x -= overlap * nx; p->y -= overlap * ny; p->z -= overlap * nz;
                
                float rvx = particles[other].vx - p->vx;
                float rvy = particles[other].vy - p->vy;
                float rvz = particles[other].vz - p->vz;
                float vn = rvx*nx + rvy*ny + rvz*nz;
                if (vn < 0) {
                    float j = -1.5f * vn;
                    p->vx -= j * nx * 0.5f; p->vy -= j * ny * 0.5f; p->vz -= j * nz * 0.5f;
                }
            }
        }
    } else {
        float dx = n->x - p->x; float dy = n->y - p->y; float dz = n->z - p->z;
        if ((dx*dx + dy*dy + dz*dz) > 400.0f && n->particle_count > 1) return;
        for(int i=0; i<8; i++) if(n->children[i] != -1) resolve_collisions(particles, p_idx, n->children[i]);
    }
}

// --- Adaptive Logic ---
void assign_timesteps(Particle* particles, SimConfig config) {
    int i;
    float eta = 0.5f; 
    
    #pragma omp parallel for private(i)
    for (i = 0; i < config.particle_count; i++) {
        float acc2 = particles[i].ax * particles[i].ax + 
                     particles[i].ay * particles[i].ay + 
                     particles[i].az * particles[i].az;
        
        float acc = sqrtf(acc2);
        if (acc < 1e-6f) acc = 1e-6f;

        float ideal_dt = eta * sqrtf(config.adaptive_err / acc);
        if (ideal_dt > config.dt) ideal_dt = config.dt;

        float ratio = config.dt / ideal_dt;
        int level = (int)ceilf(log2f(ratio));
        
        if (level < 0) level = 0;
        if (level > config.max_level) level = config.max_level;

        particles[i].level = level;
        particles[i].current_dt = config.dt / (float)(1 << level);
    }
}

// --- LEAPFROG FUNCTION ---
void leapfrog_active_step(Particle* particles, int p_idx, float dt, int root_idx, SimConfig config) {
    Particle* p = &particles[p_idx];

    // Kick 1
    p->vx += p->ax * (dt * 0.5f);
    p->vy += p->ay * (dt * 0.5f);
    p->vz += p->az * (dt * 0.5f);

    // Drift
    p->x += p->vx * dt;
    p->y += p->vy * dt;
    p->z += p->vz * dt;

    // Force Calculation
    float fx = 0, fy = 0, fz = 0;
    calculate_force(p, root_idx, config, &fx, &fy, &fz);

    p->ax = fx / p->mass;
    p->ay = fy / p->mass;
    p->az = fz / p->mass;

    // Kick 2
    p->vx += p->ax * (dt * 0.5f);
    p->vy += p->ay * (dt * 0.5f);
    p->vz += p->az * (dt * 0.5f);

    resolve_collisions(particles, p_idx, root_idx);
}