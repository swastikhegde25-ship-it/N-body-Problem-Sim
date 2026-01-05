#ifndef PHYSICS_H
#define PHYSICS_H

#include <stdint.h>

#ifdef _WIN32
    #define EXPORT __declspec(dllexport)
#else
    #define EXPORT
#endif

// Updated Struct with Morton Index
typedef struct {
    float x, y, z;
    float vx, vy, vz;
    float mass;
    uint64_t morton_index; // <-- Added for sorting
} Particle;

typedef struct {
    int count;
    float dt;
} SimConfig;

EXPORT void init_simulation();
EXPORT void step_simulation(Particle* particles, SimConfig config);
EXPORT void render_cpu(Particle* particles, int count, uint8_t* pixels, int width, int height);

#endif