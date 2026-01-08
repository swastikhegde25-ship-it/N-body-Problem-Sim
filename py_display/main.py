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
        int particle_count;
        float dt;
        float G;
        float softening;
        float world_size;
        float theta;
    } SimConfig;

    void init_simulation();
    void step_simulation(Particle* particles, SimConfig config);
    void render_cpu(Particle* particles, int count, uint8_t* pixels, int width, int height);
""")

dll_path = os.path.join(os.path.dirname(__file__), "..", "c_core", "physics.dll")
if not os.path.exists(dll_path):
    print("DLL not found. Compile c_core/physics.c")
    exit()
lib = ffi.dlopen(dll_path)

# --- Init ---
WIDTH, HEIGHT = 1000, 800
NUM_PARTICLES = 100000

particle_dtype = np.dtype([
    ('x', 'f4'), ('y', 'f4'), ('z', 'f4'),
    ('vx', 'f4'), ('vy', 'f4'), ('vz', 'f4'),
    ('mass', 'f4'),
    ('morton', 'u8')
], align=True)

particles_np = np.zeros(NUM_PARTICLES, dtype=particle_dtype)
particles_ptr = ffi.cast("Particle*", particles_np.ctypes.data)

def create_sphere(start, count, cx, cy, cz, r, mass, vx, vy, vz):
    rad = np.random.uniform(0, 1, count) ** (1/3) * r
    theta = np.random.uniform(0, 2*np.pi, count)
    phi = np.random.uniform(0, np.pi, count)
    x = rad * np.sin(phi) * np.cos(theta)
    y = rad * np.sin(phi) * np.sin(theta)
    z = rad * np.cos(phi)
    
    end = start + count
    particles_np['x'][start:end] = x + cx
    particles_np['y'][start:end] = y + cy
    particles_np['z'][start:end] = z + cz
    particles_np['vx'][start:end] = vx
    particles_np['vy'][start:end] = vy
    particles_np['vz'][start:end] = vz
    particles_np['mass'][start:end] = mass / count

print("Generating 2 Planets...")
create_sphere(0, NUM_PARTICLES//2, -400, 0, 0, 200, 100000, 15, 10, 0)
create_sphere(NUM_PARTICLES//2, NUM_PARTICLES//2, 400, 0, 0, 200, 100000, -15, -10, 0)

config = ffi.new("SimConfig*")
config.particle_count = NUM_PARTICLES
config.dt = 0.05
config.G = 1.0
config.softening = 2.0
config.world_size = 4000.0
config.theta = 0.5

lib.init_simulation()

# --- Loop ---
pygame.init()
screen = pygame.display.set_mode((WIDTH, HEIGHT))
clock = pygame.time.Clock()

pixel_buffer = np.zeros((WIDTH * HEIGHT * 3), dtype=np.uint8)
pixel_ptr = ffi.cast("uint8_t*", pixel_buffer.ctypes.data)

running = True
while running:
    lib.step_simulation(particles_ptr, config[0])
    lib.render_cpu(particles_ptr, NUM_PARTICLES, pixel_ptr, WIDTH, HEIGHT)
    
    surface_array = pixel_buffer.reshape((HEIGHT, WIDTH, 3))
    surface = pygame.surfarray.make_surface(surface_array.swapaxes(0, 1))
    screen.blit(surface, (0, 0))
    
    pygame.display.set_caption(f"Commit 3: CPU Physics (Serial) | FPS: {clock.get_fps():.2f}")
    pygame.display.flip()
    clock.tick(60)

    for event in pygame.event.get():
        if event.type == pygame.QUIT: running = False

pygame.quit()