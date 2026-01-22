import pygame
import numpy as np
import cffi
import os
import math
import matplotlib.pyplot as plt
from collections import deque

# used github copilot inline code suggestion to speed up development 
# but all the code logic, structure, and implementation is my own work.

# --- Setup ---
ffi = cffi.FFI()
ffi.cdef("""
    typedef struct {
        float x, y, z;
        float vx, vy, vz;
        float ax, ay, az;
        float mass;
        int level;
        float current_dt;
        uint64_t morton_index;
    } Particle;

    typedef struct {
        int particle_count;
        float dt;
        float G;
        float softening;
        float world_size;
        float theta;
        int max_level;
        float adaptive_err;
    } SimConfig;
         
    typedef struct {
        double total_energy;
        double kinetic;
        double potential;
    } EnergyStats;

    void init_simulation();
    void step_simulation(Particle* particles, SimConfig config);
    void get_energy_stats(Particle* particles, SimConfig config, EnergyStats* stats);
    
    // UPDATED: Now takes a matrix pointer
    void render_cpu(Particle* particles, int count, uint8_t* pixels, 
                    int width, int height, 
                    float* rot_matrix, float zoom_factor);
""")

dll_path = os.path.join(os.path.dirname(__file__), "..", "c_core", "physics.dll")
if not os.path.exists(dll_path):
    print("DLL not found. Build it!")
    exit()
lib = ffi.dlopen(dll_path)

WIDTH, HEIGHT = 1000, 800
NUM_PARTICLES = 100000

particle_dtype = np.dtype([
    ('x', 'f4'), ('y', 'f4'), ('z', 'f4'),
    ('vx', 'f4'), ('vy', 'f4'), ('vz', 'f4'),
    ('ax', 'f4'), ('ay', 'f4'), ('az', 'f4'),
    ('mass', 'f4'),
    ('level', 'i4'),
    ('current_dt', 'f4'),
    ('morton', 'u8')
], align=True)

particles_np = np.zeros(NUM_PARTICLES, dtype=particle_dtype)
particles_ptr = ffi.cast("Particle*", particles_np.ctypes.data)

# Init Planets
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
    # Initialize acc to 0
    particles_np['ax'][start:end] = 0
    particles_np['ay'][start:end] = 0
    particles_np['az'][start:end] = 0

print("Generating Planets...")
create_sphere(0, NUM_PARTICLES, 0, 0, 0, 200, 100000, 0, 0, 0)
# create_sphere(0, NUM_PARTICLES//2, -400, 0, 0, 200, 100000, 15, 10, 0)
# create_sphere(NUM_PARTICLES//2, NUM_PARTICLES//2, 400, 0, 0, 200, 100000, -15, -10, 0)

config = ffi.new("SimConfig*")
config.particle_count = NUM_PARTICLES
config.dt = 0.1
config.G = 1.0
config.softening = 2.0
config.world_size = 4000.0
config.theta = 0.5
# Block Adaptive Config
config.max_level = 4       # Max hierarchy depth 
config.adaptive_err = 10.0 # Control factor for dt calculation

lib.init_simulation()

pygame.init()
screen = pygame.display.set_mode((WIDTH, HEIGHT))
clock = pygame.time.Clock()
pixel_buffer = np.zeros((WIDTH * HEIGHT * 3), dtype=np.uint8)
pixel_ptr = ffi.cast("uint8_t*", pixel_buffer.ctypes.data)

# --------------- Graph Setup ----------------
plt.ion()
fig, ax = plt.subplots(figsize=(5, 3))
history_len = 200
t_data = deque(maxlen=history_len)
ke_data = deque(maxlen=history_len)
pe_data = deque(maxlen=history_len)
te_data = deque(maxlen=history_len)

line_ke, = ax.plot([], [], label='KE', color='r')
line_pe, = ax.plot([], [], label='PE', color='b')
line_te, = ax.plot([], [], label='Total', color='g')
ax.legend()
ax.set_title("System Energy")
ax.grid(True, alpha=0.3)

energy_stats = ffi.new("EnergyStats*")
frame_count = 0

# --- ROTATION MATRIX LOGIC ---
# Identity Matrix (No rotation)
# [ R00 R01 R02 ]
# [ R10 R11 R12 ]
# [ R20 R21 R22 ]
rotation_matrix = np.eye(3, dtype=np.float32)

def rotate_x(angle):
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]], dtype=np.float32)

def rotate_y(angle):
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]], dtype=np.float32)

cam_zoom = 300.0
last_mouse_pos = (0, 0)
dragging_left = False
dragging_right = False
paused = False

print("Controls:")
print("  [Left Drag]   Pitch (Vertical Axis Rotation)")
print("  [Right Drag]  Yaw   (Horizontal Axis Rotation)")
print("  [Scroll]      Zoom")
print("  [Space]       Pause")
print("  [R]           Reset")

running = True
while running:
    if not paused:
        lib.step_simulation(particles_ptr, config[0])
        
        # Update Graph every 10 frames to save FPS
        if frame_count % 10 == 0:
            lib.get_energy_stats(particles_ptr, config[0], energy_stats)
            t_data.append(frame_count)
            ke_data.append(energy_stats.kinetic)
            pe_data.append(energy_stats.potential)
            te_data.append(energy_stats.total_energy)
            
            line_ke.set_data(t_data, ke_data)
            line_pe.set_data(t_data, pe_data)
            line_te.set_data(t_data, te_data)
            
            ax.relim()
            ax.autoscale_view()
            plt.pause(0.001)

    frame_count += 1

    # Pass the matrix to C
    # flatten() converts 3x3 to 9 floats
    rot_ptr = ffi.cast("float*", rotation_matrix.ctypes.data)
    
    lib.render_cpu(particles_ptr, NUM_PARTICLES, pixel_ptr, WIDTH, HEIGHT, 
                   rot_ptr, cam_zoom)
    
    surface_array = pixel_buffer.reshape((HEIGHT, WIDTH, 3))
    surface = pygame.surfarray.make_surface(surface_array.swapaxes(0, 1))
    screen.blit(surface, (0, 0))
    
    state = "PAUSED" if paused else "RUNNING"
    # adding a fps counter to the window title
    fps = clock.get_fps()
    pygame.display.set_caption(f"Sim | {state} | Particles: {NUM_PARTICLES} | FPS: {fps:.2f}")

    pygame.display.flip()
    clock.tick(60)

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
        elif event.type == pygame.KEYDOWN:
            if event.key == pygame.K_ESCAPE: running = False
            if event.key == pygame.K_SPACE: paused = not paused
            if event.key == pygame.K_r: 
                rotation_matrix = np.eye(3, dtype=np.float32)
                cam_zoom = 300.0
                
        elif event.type == pygame.MOUSEBUTTONDOWN:
            if event.button == 1: # Left Click
                dragging_left = True
                last_mouse_pos = event.pos
            elif event.button == 3: # Right Click
                dragging_right = True
                last_mouse_pos = event.pos
            elif event.button == 4: cam_zoom *= 1.1
            elif event.button == 5: cam_zoom /= 1.1

        elif event.type == pygame.MOUSEBUTTONUP:
            if event.button == 1: dragging_left = False
            if event.button == 3: dragging_right = False

        elif event.type == pygame.MOUSEMOTION:
            dx = event.pos[0] - last_mouse_pos[0]
            dy = event.pos[1] - last_mouse_pos[1]
            last_mouse_pos = event.pos
            
            sensitivity = 0.005
            
            # Left Click: Vertical Mouse Movement -> Rotate around X axis (Pitch)
            if dragging_left:
                rot = rotate_x(dy * sensitivity)
                # Apply rotation relative to current view
                rotation_matrix = np.dot(rot, rotation_matrix)
                
            # Right Click: Horizontal Mouse Movement -> Rotate around Y axis (Yaw)
            if dragging_right:
                rot = rotate_y(dx * sensitivity)
                rotation_matrix = np.dot(rot, rotation_matrix)
