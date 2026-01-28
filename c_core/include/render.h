#pragma once
#include "common.h"

EXPORT void render_cpu(Particle* particles, int count, uint8_t* pixels, 
                       int width, int height, 
                       float* rot_matrix, float zoom_factor, 
                       float tx, float ty, float tz);