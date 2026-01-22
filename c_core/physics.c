#include "physics.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h> 

/* used github copilot inline code suggestion to speed up development 
but all the code logic, structure, and implementation is my own work.
Took help for some parts of the logic from stackoverflow */

#undef MAX_NODES
#define MAX_NODES 2000000

Node node_pool[MAX_NODES];
int node_count = 0;

// --- Morton & Sort (Standard) ---
uint64_t expand_bits(uint32_t v) {
    uint64_t x = v & 0x1fffff;
    x = (x | x << 32) & 0x1f00000000ffff;
    x = (x | x << 16) & 0x1f0000ff0000ff;
    x = (x | x << 8)  & 0x100f00f00f00f00f;
    x = (x | x << 4)  & 0x10c30c30c30c30c3;
    x = (x | x << 2)  & 0x1249249249249249;
    return x;
}

uint64_t morton_3d(float x, float y, float z, float min_c, float max_c) {
    float range = max_c - min_c;
    float nx = (x - min_c) / range;
    float ny = (y - min_c) / range;
    float nz = (z - min_c) / range;
    if(nx < 0) nx = 0; if(nx > 1) nx = 1;
    if(ny < 0) ny = 0; if(ny > 1) ny = 1;
    if(nz < 0) nz = 0; if(nz > 1) nz = 1;
    uint32_t ix = (uint32_t)(nx * 2097151.0f);
    uint32_t iy = (uint32_t)(ny * 2097151.0f);
    uint32_t iz = (uint32_t)(nz * 2097151.0f);
    return (expand_bits(ix) << 2) | (expand_bits(iy) << 1) | expand_bits(iz);
}

int compare_particles(const void* a, const void* b) {
    const Particle* pa = (const Particle*)a;
    const Particle* pb = (const Particle*)b;
    if (pa->morton_index < pb->morton_index) return -1;
    if (pa->morton_index > pb->morton_index) return 1;
    return 0;
}

// --- Tree & Physics (Standard) ---
void reset_tree() { node_count = 0; }

int alloc_node() {
    int idx = node_count++;
    if (idx >= MAX_NODES) return -1;
    Node* n = &node_pool[idx];
    n->mass = 0; n->x = 0; n->y = 0; n->z = 0;
    n->first_particle = -1; n->particle_count = 0;
    for(int i=0; i<8; i++) n->children[i] = -1;
    return idx;
}

void insert_node(int node_idx, int p_idx, Particle* p) {
    Node* n = &node_pool[node_idx];
    n->mass += p->mass;
    n->x += p->x * p->mass; n->y += p->y * p->mass; n->z += p->z * p->mass;
    n->particle_count++;
    if (n->particle_count == 1) n->first_particle = p_idx;

    float mid_x = (n->min_x + n->max_x) * 0.5f;
    float mid_y = (n->min_y + n->max_y) * 0.5f;
    float mid_z = (n->min_z + n->max_z) * 0.5f;
    int octant = 0;
    if (p->x > mid_x) octant |= 1;
    if (p->y > mid_y) octant |= 2;
    if (p->z > mid_z) octant |= 4;

    if (n->children[octant] == -1) {
        int child_idx = alloc_node();
        if(child_idx == -1) return; // Pool full, child not created. Mass stays in parent.
        Node* child = &node_pool[child_idx];
        child->min_x = (octant & 1) ? mid_x : n->min_x;
        child->max_x = (octant & 1) ? n->max_x : mid_x;
        child->min_y = (octant & 2) ? mid_y : n->min_y;
        child->max_y = (octant & 2) ? n->max_y : mid_y;
        child->min_z = (octant & 4) ? mid_z : n->min_z;
        child->max_z = (octant & 4) ? n->max_z : mid_z;
        n->children[octant] = child_idx;
    }
    // Limit depth to prevent stack overflow/infinite recursion on exact overlaps
    if ((n->max_x - n->min_x) < 0.5f) return;
    insert_node(n->children[octant], p_idx, p);
}

void finalize_nodes(int node_idx) {
    Node* n = &node_pool[node_idx];
    if (n->mass > 0) {
        n->x /= n->mass; n->y /= n->mass; n->z /= n->mass;
    }
    for(int i=0; i<8; i++) if(n->children[i] != -1) finalize_nodes(n->children[i]);
}

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

    int root = alloc_node();
    node_pool[root].min_x = -config.world_size; node_pool[root].max_x = config.world_size;
    node_pool[root].min_y = -config.world_size; node_pool[root].max_y = config.world_size;
    node_pool[root].min_z = -config.world_size; node_pool[root].max_z = config.world_size;

    for (i = 0; i < config.particle_count; i++) insert_node(root, i, &particles[i]);
    finalize_nodes(root);

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
                leapfrog_active_step(particles, i, particles[i].current_dt, root, config);
            }
        }
    }
}

// --- MATRIX RENDERER (CAD Like) ---
EXPORT void render_cpu(Particle* particles, int count, uint8_t* pixels, 
                       int width, int height, 
                       float* rot_matrix, float zoom_factor) {
    memset(pixels, 0, width * height * 3);

    float r00 = rot_matrix[0]; float r01 = rot_matrix[1]; float r02 = rot_matrix[2];
    float r10 = rot_matrix[3]; float r11 = rot_matrix[4]; float r12 = rot_matrix[5];
    float r20 = rot_matrix[6]; float r21 = rot_matrix[7]; float r22 = rot_matrix[8];

    float brightness_scale = zoom_factor / 300.0f;
    if (brightness_scale < 0.5f) brightness_scale = 0.5f;
    int add_b = (int)(5 * brightness_scale);
    int add_g = (int)(10 * brightness_scale);
    int add_r = (int)(20 * brightness_scale);
    if (add_b > 255) add_b = 255; if (add_g > 255) add_g = 255; if (add_r > 255) add_r = 255;

    int i;
    #pragma omp parallel for private(i)
    for (i = 0; i < count; i++) {
        float x = particles[i].x;
        float y = particles[i].y;
        float z = particles[i].z;

        float rx = x * r00 + y * r01 + z * r02;
        float ry = x * r10 + y * r11 + z * r12;
        float rz = x * r20 + y * r21 + z * r22;

        float dist = 1000.0f; 
        float final_z = dist + rz;

        if (final_z < 10.0f) continue;

        float scale = zoom_factor / final_z;
        int sx = (int)(rx * scale + (width / 2));
        int sy = (int)(ry * scale + (height / 2));

        if (sx >= 0 && sx < width && sy >= 0 && sy < height) {
            int index = (sy * width + sx) * 3;
            #pragma omp atomic
            pixels[index] += add_b;
            #pragma omp atomic
            pixels[index+1] += add_g;
            #pragma omp atomic
            pixels[index+2] += add_r;
        }
    }
}