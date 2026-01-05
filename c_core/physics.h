#ifndef PHYSICS_H
#define PHYSICS_H

#include <stdint.h>

#ifdef _WIN32
    #define EXPORT __declspec(dllexport)
#else
    #define EXPORT
#endif

// Simple Particle Struct (No padding needed yet)
typedef struct {
    float x, y, z;
    float vx, vy, vz;
    float mass;
} Particle;

// Config Struct
typedef struct {
    int count;
    float dt;
} SimConfig;

// Exported Functions
EXPORT void init_simulation();
EXPORT void step_simulation(Particle* particles, SimConfig config);
EXPORT void render_cpu(Particle* particles, int count, uint8_t* pixels, int width, int height);

#endif