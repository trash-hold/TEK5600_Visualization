#pragma once

#include <random>

#include "vec2.h"
#include "data_structures.h"
#include "integrators.h"
#include "helpers.h"

PixelPlane simpleLIC(const RawData* data, const float& step_size ,const size_t& kernel_length);