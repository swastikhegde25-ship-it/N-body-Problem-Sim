#include "morton.h"

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