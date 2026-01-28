#pragma once
#include "common.h"

uint64_t expand_bits(uint32_t v);
uint64_t morton_3d(float x, float y, float z, float min_c, float max_c);
int compare_particles(const void* a, const void* b);