#pragma once

#include <vector>
#include <deque>

#include "vec2.h"
#include "data_structures.h"

struct Line{
    std::vector<Particle> line;
};

Particle euler_step(const Vec2& pos, const RawData* data, const float& step_size);

std::vector<Line> eulerIntegrator(const std::vector<Particle>* seeds, const RawData* data, const float &step_size, const size_t &max_steps);

Particle rk4_step(const Vec2& pos, const RawData* data, const float& step_size);

std::vector<Line> rk4Integrator(const std::vector<Particle>* seeds, const RawData* data, const float &step_size, const size_t &max_steps);