#ifndef PHYSICS_H
#define PHYSICS_H

#include <stdint.h>

#ifdef _WIN32
    #define EXPORT __declspec(dllexport)
#else
    #define EXPORT
#endif

#define MAX_NODES 400000 

typedef struct {
    float x, y, z;
    float vx, vy, vz;
    float mass;
    uint64_t morton_index;
} Particle;

// Octree Node
typedef struct {
    float mass;
    float x, y, z;          // Center of Mass
    float min_x, max_x;     // Bounds
    float min_y, max_y;
    float min_z, max_z;
    int children[8];
    int first_particle;     // For collision tracking
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

EXPORT void init_simulation();
EXPORT void step_simulation(Particle* particles, SimConfig config);
EXPORT void render_cpu(Particle* particles, int count, uint8_t* pixels, int width, int height);

#endif