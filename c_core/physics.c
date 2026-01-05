#include "physics.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// --- Morton Logic ---
// Expands a 21-bit integer into 64 bits by inserting 2 zeros after each bit.
uint64_t expand_bits(uint32_t v) {
    uint64_t x = v & 0x1fffff;
    x = (x | x << 32) & 0x1f00000000ffff;
    x = (x | x << 16) & 0x1f0000ff0000ff;
    x = (x | x << 8)  & 0x100f00f00f00f00f;
    x = (x | x << 4)  & 0x10c30c30c30c30c3;
    x = (x | x << 2)  & 0x1249249249249249;
    return x;
}

// Calculates Z-Order curve index for 3D coordinates
uint64_t morton_3d(float x, float y, float z) {
    // 1. Normalize to World Bounds (Assume -1000 to 1000)
    float offset = 1000.0f;
    float range = 2000.0f;
    
    float ox = x + offset;
    float oy = y + offset;
    float oz = z + offset;

    // 2. Clamp
    if(ox < 0) ox = 0; if(ox > range) ox = range;
    if(oy < 0) oy = 0; if(oy > range) oy = range;
    if(oz < 0) oz = 0; if(oz > range) oz = range;

    // 3. Map to Integer Grid (21 bits per axis)
    uint32_t ix = (uint32_t)((ox / range) * 2097151.0f);
    uint32_t iy = (uint32_t)((oy / range) * 2097151.0f);
    uint32_t iz = (uint32_t)((oz / range) * 2097151.0f);

    // 4. Interleave bits
    return (expand_bits(ix) << 2) | (expand_bits(iy) << 1) | expand_bits(iz);
}

// Comparator for QuickSort
int compare_particles(const void* a, const void* b) {
    const Particle* pa = (const Particle*)a;
    const Particle* pb = (const Particle*)b;
    if (pa->morton_index < pb->morton_index) return -1;
    if (pa->morton_index > pb->morton_index) return 1;
    return 0;
}

// --- Main Functions ---

EXPORT void init_simulation() {
    printf("[C] Physics Initialized. Z-Sorting Enabled.\n");
}

EXPORT void step_simulation(Particle* particles, SimConfig config) {
    // 1. Update Positions (Simple Rotation)
    float cos_t = cosf(config.dt);
    float sin_t = sinf(config.dt);

    for (int i = 0; i < config.count; i++) {
        float nx = particles[i].x * cos_t - particles[i].z * sin_t;
        float nz = particles[i].x * sin_t + particles[i].z * cos_t;
        particles[i].x = nx;
        particles[i].z = nz;

        // 2. Calculate Morton Code
        particles[i].morton_index = morton_3d(particles[i].x, particles[i].y, particles[i].z);
    }

    // 3. SORT the Array physically in RAM
    qsort(particles, config.count, sizeof(Particle), compare_particles);
}

void project_particle(float x, float y, float z, int* sx, int* sy, int sw, int sh) {
    float fov = 300.0f;
    float viewer_dist = 1000.0f;
    float scale = fov / (viewer_dist + z);
    *sx = (int)(x * scale + (sw / 2));
    *sy = (int)(y * scale + (sh / 2));
}

EXPORT void render_cpu(Particle* particles, int count, uint8_t* pixels, int width, int height) {
    memset(pixels, 0, width * height * 3);

    for (int i = 0; i < count; i++) {
        int sx, sy;
        project_particle(particles[i].x, particles[i].y, particles[i].z, &sx, &sy, width, height);
        
        if (sx >= 0 && sx < width && sy >= 0 && sy < height) {
            int index = (sy * width + sx) * 3;
            
            // --- DEBUG VISUALIZATION ---
            // Color based on "Rank" in the sorted array.
            // If sorted correctly, 'i' correlates to 3D position.
            // This creates the Rainbow Gradient effect.
            
            float t = (float)i / (float)count; // 0.0 to 1.0
            
            uint8_t r = (uint8_t)(t * 255.0f);
            uint8_t b = (uint8_t)((1.0f - t) * 255.0f);
            uint8_t g = (uint8_t)(sinf(t * 3.14159f) * 255.0f);

            pixels[index]   = b;
            pixels[index+1] = g;
            pixels[index+2] = r;
        }
    }
}