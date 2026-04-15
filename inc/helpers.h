#pragma once

#include <iostream>
#include <H5Cpp.h>

#include "vec2.h"
#include "data_structures.h"

// Function to read the HDF5 file -- there is assumed structure of the file and equality of dimensions for both x and y components
void readH5File(const std::string& filename, RawData* data);

float bi_interpolate(float v00, float v01, float v10, float v11, float x_frac, float y_frac);

Vec2 bi_interpolate(const RawData* data, const Vec2& pos);

Vec2 bi_interpolate(const PixelPlane* data, const Vec2& pos);