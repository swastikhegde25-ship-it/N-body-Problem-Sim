#include "physics.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

EXPORT void init_simulation() {
    printf("[C] Simulation Initialized.\n");
}

EXPORT void step_simulation(Particle* particles, SimConfig config) {
    // Simple Rotation Logic (No Gravity yet)
    // This proves the C code is modifying memory
    float cos_t = cosf(config.dt);
    float sin_t = sinf(config.dt);

    for (int i = 0; i < config.count; i++) {
        float nx = particles[i].x * cos_t - particles[i].z * sin_t;
        float nz = particles[i].x * sin_t + particles[i].z * cos_t;
        particles[i].x = nx;
        particles[i].z = nz;
    }
}

// Helper to project 3D to 2D
void project_particle(float x, float y, float z, int* sx, int* sy, int sw, int sh) {
    float fov = 300.0f;
    float viewer_dist = 1000.0f;
    float scale = fov / (viewer_dist + z);
    *sx = (int)(x * scale + (sw / 2));
    *sy = (int)(y * scale + (sh / 2));
}

EXPORT void render_cpu(Particle* particles, int count, uint8_t* pixels, int width, int height) {
    // Clear screen to black
    memset(pixels, 0, width * height * 3);

    for (int i = 0; i < count; i++) {
        int sx, sy;
        project_particle(particles[i].x, particles[i].y, particles[i].z, &sx, &sy, width, height);
        
        // Draw simple white pixel if on screen
        if (sx >= 0 && sx < width && sy >= 0 && sy < height) {
            int index = (sy * width + sx) * 3;
            pixels[index] = 255;   // B
            pixels[index+1] = 255; // G
            pixels[index+2] = 255; // R
        }
    }
}