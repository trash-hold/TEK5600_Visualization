#pragma once

#include <iostream>
#include <H5Cpp.h>

#include "vec2.h"
#include "data_structures.h"

// Function to read the HDF5 file -- there is assumed structure of the file and equality of dimensions for both x and y components
void readH5File(const std::string& filename, RawData* data);