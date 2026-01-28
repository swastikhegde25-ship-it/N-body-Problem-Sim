#pragma once
#include "common.h"
#include "morton.h"
#include "tree.h"
#include "engine.h"
#include "energy.h"
#include "render.h"

EXPORT void init_simulation();
EXPORT void step_simulation(Particle* particles, SimConfig config);