#pragma once

#include <stdint.h>

#ifdef _WIN32
    #define EXPORT __declspec(dllexport)
#else
    #define EXPORT
#endif

#define MAX_NODES 400000 

// --- Structures ---
typedef struct {
    float x, y, z;
    float vx, vy, vz;
    float mass;
    float padding; // Keeping padding for alignment safety
    uint64_t morton_index;
} Particle;

typedef struct {
    float mass;
    float x, y, z;
    float min_x, max_x;
    float min_y, max_y;
    float min_z, max_z;
    int children[8];
    int first_particle;
    int particle_count;
} Node;

typedef struct {
    int particle_count;
    float dt;
    float G;
    float softening;
    float world_size;
    float theta;
} SimConfig;

// --- Functions ---
EXPORT void init_simulation();
EXPORT void step_simulation(Particle* particles, SimConfig config);

// accepts a 3x3 Rotation Matrix (array of 9 floats)
EXPORT void render_cpu(Particle* particles, int count, uint8_t* pixels, 
                       int width, int height, 
                       float* rot_matrix, float zoom_factor);