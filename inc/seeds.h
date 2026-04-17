#pragma once

#include <vector>
#include <random>
#include <deque>
#include <iostream>

#include "vec2.h"
#include "data_structures.h"
#include "integrators.h"
#include "helpers.h"

std::vector<Particle> getRandomSeed(RawData* const data, float particle_percentage);

std::vector<Particle> getUniformSeed(RawData* const data, size_t cell_distance);

bool growFieldLine(SpatialHash& hash, const RawData* data, std::vector<Line>& field_lines, const Vec2& seed, const uint& lineID, const float &step_size, const size_t &max_steps);

std::vector<Line> getEvenSeed(RawData* const data, const float& line_distance, const float& step_size, const size_t& max_steps);