#pragma once

#include <random>

#include "vec2.h"
#include "data_structures.h"
#include "integrators.h"
#include "helpers.h"

template <typename F>
PixelPlane simpleLIC(const RawData* data, const float& step_size, const size_t& kernel_length, F integrator);
template <typename F2>
PixelPlane simpleScaledLIC(const RawData* data, const size_t &scale_factor, const float& step_size, const size_t& kernel_length, F2 integrator);