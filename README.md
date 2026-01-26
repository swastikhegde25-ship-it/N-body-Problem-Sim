# N-Body Gravity Simulation

A high-performance gravity simulation capable of rendering **100k + particles** in real-time. This project uses a **Hybrid Architecture**:
C is used for heavy physics calculations, and Python is used for visualization.

##  Features & Technical Details

*   **Hybrid Engine:** Python controls the simulation, C calculates the math.
*   **Barnes-Hut Algorithm:** Optimized physics using an Octree (N log N complexity).    
*   **Parallel Processing:** Uses OpenMP to run physics on all CPU cores.                
*   **Hierarchy Adaptive Timestep:** "dt" changes according to particle acceleration.                  
*   **Interactive:** 360-degree camera control with Zoom and Rotation.

##  Controls

*   **Spacebar:** Pause / Resume Simulation
*   **Left Click + Drag:** Rotate Camera (Pitch) & Rotate Camera (Yaw)
*   **Scroll Wheel:**  Zoom In / Out  
*   **Scroll Wheel Click:** Pan (Drag on plane) 
*   **R Key:** Reset Camera
*   **ESC:** Quit

##  How to Run

1.  **Install Requirements:**
    ```bash
    pip install pygame numpy cffi moderngl
    ```

2.  **Compile the Physics Engine:**
    *   In VS Code, press `Ctrl + Shift + B` to build the C DLL.

3.  **Run the Simulation:**
    ```bash
    python main.py
    ```
