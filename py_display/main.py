import pygame
import numpy as np
import cffi
import os

# --- CFFI Setup ---
ffi = cffi.FFI()
ffi.cdef("""
    typedef struct {
        float x, y, z;
        float vx, vy, vz;
        float mass;
        uint64_t morton_index;
    } Particle;

    typedef struct {
        int count;
        float dt;
    } SimConfig;

    void init_simulation();
    void step_simulation(Particle* particles, SimConfig config);
    void render_cpu(Particle* particles, int count, uint8_t* pixels, int width, int height);
""")

# Path correction for your folder structure
dll_path = os.path.join(os.path.dirname(__file__), "..", "c_core", "physics.dll")
if not os.path.exists(dll_path):
    print(f"DLL not found at: {dll_path}")
    print("Please compile physics.c in the c_core folder!")
    exit()
lib = ffi.dlopen(dll_path)

# --- Init Data ---
WIDTH, HEIGHT = 1000, 800
NUM_PARTICLES = 100000

# Updated Dtype with Morton (u8 = uint64)
particle_dtype = np.dtype([
    ('x', 'f4'), ('y', 'f4'), ('z', 'f4'),
    ('vx', 'f4'), ('vy', 'f4'), ('vz', 'f4'),
    ('mass', 'f4'),
    ('morton', 'u8') 
], align=True)

particles_np = np.zeros(NUM_PARTICLES, dtype=particle_dtype)
particles_ptr = ffi.cast("Particle*", particles_np.ctypes.data)

# Random Sphere Init
r = np.random.uniform(0, 200, NUM_PARTICLES)
theta = np.random.uniform(0, 2*np.pi, NUM_PARTICLES)
phi = np.random.uniform(0, np.pi, NUM_PARTICLES)

particles_np['x'] = r * np.sin(phi) * np.cos(theta)
particles_np['y'] = r * np.sin(phi) * np.sin(theta)
particles_np['z'] = r * np.cos(phi)

# Config
config = ffi.new("SimConfig*")
config.count = NUM_PARTICLES
config.dt = 0.01

lib.init_simulation()

# --- Loop ---
pygame.init()
screen = pygame.display.set_mode((WIDTH, HEIGHT))
clock = pygame.time.Clock()

pixel_buffer = np.zeros((WIDTH * HEIGHT * 3), dtype=np.uint8)
pixel_ptr = ffi.cast("uint8_t*", pixel_buffer.ctypes.data)

running = True
while running:
    # 1. Physics (Sort + Spin)
    lib.step_simulation(particles_ptr, config[0])
    
    # 2. Render
    lib.render_cpu(particles_ptr, NUM_PARTICLES, pixel_ptr, WIDTH, HEIGHT)
    
    # 3. Display
    surface_array = pixel_buffer.reshape((HEIGHT, WIDTH, 3))
    surface = pygame.surfarray.make_surface(surface_array.swapaxes(0, 1))
    screen.blit(surface, (0, 0))
    
    pygame.display.set_caption(f"Commit 2: Z-Morton Sorting | FPS: {clock.get_fps():.2f}")
    pygame.display.flip()
    clock.tick(60)

    for event in pygame.event.get():
        if event.type == pygame.QUIT: running = False

pygame.quit()