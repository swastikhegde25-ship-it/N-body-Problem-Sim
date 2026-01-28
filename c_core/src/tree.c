#include "tree.h"
#include <stdlib.h>

Node node_pool[MAX_NODES];
int node_count = 0;
int global_root = -1;

// --- Tree & Physics (Standard) ---
void reset_tree() { node_count = 0; }

int alloc_node() {
    int idx = node_count++;
    if (idx >= MAX_NODES) return -1;
    Node* n = &node_pool[idx];
    n->mass = 0; n->x = 0; n->y = 0; n->z = 0;
    n->first_particle = -1; n->particle_count = 0;
    for(int i=0; i<8; i++) n->children[i] = -1;
    return idx;
}

void insert_node(int node_idx, int p_idx, Particle* p) {
    Node* n = &node_pool[node_idx];
    n->mass += p->mass;
    n->x += p->x * p->mass; n->y += p->y * p->mass; n->z += p->z * p->mass;
    n->particle_count++;
    if (n->particle_count == 1) n->first_particle = p_idx;

    float mid_x = (n->min_x + n->max_x) * 0.5f;
    float mid_y = (n->min_y + n->max_y) * 0.5f;
    float mid_z = (n->min_z + n->max_z) * 0.5f;
    int octant = 0;
    if (p->x > mid_x) octant |= 1;
    if (p->y > mid_y) octant |= 2;
    if (p->z > mid_z) octant |= 4;

    if (n->children[octant] == -1) {
        int child_idx = alloc_node();
        if(child_idx == -1) return; // Pool full, child not created. Mass stays in parent.
        Node* child = &node_pool[child_idx];
        child->min_x = (octant & 1) ? mid_x : n->min_x;
        child->max_x = (octant & 1) ? n->max_x : mid_x;
        child->min_y = (octant & 2) ? mid_y : n->min_y;
        child->max_y = (octant & 2) ? n->max_y : mid_y;
        child->min_z = (octant & 4) ? mid_z : n->min_z;
        child->max_z = (octant & 4) ? n->max_z : mid_z;
        n->children[octant] = child_idx;
    }
    // Limit depth to prevent stack overflow/infinite recursion on exact overlaps
    if ((n->max_x - n->min_x) < 0.5f) return;
    insert_node(n->children[octant], p_idx, p);
}

void finalize_nodes(int node_idx) {
    Node* n = &node_pool[node_idx];
    if (n->mass > 0) {
        n->x /= n->mass; n->y /= n->mass; n->z /= n->mass;
    }
    for(int i=0; i<8; i++) if(n->children[i] != -1) finalize_nodes(n->children[i]);
}