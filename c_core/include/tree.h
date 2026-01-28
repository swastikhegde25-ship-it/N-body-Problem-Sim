#pragma once
#include "common.h"

extern Node node_pool[MAX_NODES];
extern int node_count;
extern int global_root;

void reset_tree();
int alloc_node();
void insert_node(int node_idx, int p_idx, Particle* p);
void finalize_nodes(int node_idx);