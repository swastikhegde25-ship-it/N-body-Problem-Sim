#include "render.h"
#include <string.h>
#include <omp.h>

// --- MATRIX RENDERER (CAD Like) ---
EXPORT void render_cpu(Particle* particles, int count, uint8_t* pixels, 
                       int width, int height, 
                       float* rot_matrix, float zoom_factor, 
                       float tx, float ty, float tz) {
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
        // Shift particle coordinate by Target (Pivot) position first
        float dx = particles[i].x - tx;
        float dy = particles[i].y - ty;
        float dz = particles[i].z - tz;

        // Apply Rotation to the shifted coordinates
        float rx = dx * r00 + dy * r01 + dz * r02;
        float ry = dx * r10 + dy * r11 + dz * r12;
        float rz = dx * r20 + dy * r21 + dz * r22;

        float dist = 1000.0f; 
        float final_z = dist + rz;

        if (final_z < 10.0f) continue;

        float scale = zoom_factor / final_z;
        // Project. Since we already subtracted target, (0,0,0) is center screen.
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