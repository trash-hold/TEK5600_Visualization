#pragma once

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdint>
#include <stdexcept>

#include "vec2.h"

// A struct to hold the raw data read from the files
struct RawData{
    size_t n_rows, n_cols; // dimensions of the data;
    std::vector<float> x;
    std::vector<float> y;

    bool isValid(const size_t& i, const size_t& j) const {
        if (i >= n_rows || j >= n_cols || i < 0 || j < 0) {
            return false;
        }

        return true;
    }

    Vec2 at(const size_t& i, const size_t& j) const {
        if (isValid(i, j) == false) {
            throw std::out_of_range("Index out of bounds");
        }
        size_t index = i*n_cols + j;
        return Vec2(x[index], y[index]);
    }
};

struct PixelPlane{
    size_t n_rows, n_cols;
    std::vector<float> value;   // Flattened 2D array
    
    const float eps = 1e-4f;
    
    float at(size_t i, size_t j) const {
        if (i >= n_rows || j >= n_cols || i < 0 || j < 0) {
            throw std::out_of_range("Index out of bounds");
        }
        return value[i * n_cols + j];
    }

    void set(size_t i, size_t j, float val) {
        if (i >= n_rows || j >= n_cols || i < 0 || j < 0) {
            throw std::out_of_range("Index out of bounds");
        }
        value[i * n_cols + j] = val;
    }

    float bi_interpolate(const float &v00, const float& v01, const float& v10, const float& v11, const float &x_frac, const float &y_frac) const {
        // The points are interprated as follows:
        // (v10)---------(v11)
        //    |             |
        //    |        P    |
        //    |         .   |     _
        //    |             |     |
        //    |             |     |  y_frac
        // (v00)---------(v01)    -
        //
        //    |---------|  x_frac
        // 
        // And the billinear interpolation is done according to the formula found in the paper. 
        
        float v0 = v00*(1-x_frac) + v01*x_frac;
        float v1 = v10*(1-x_frac) + v11*x_frac;

        return v0*(1-y_frac) + v1*y_frac;
    }

    float interpolate(const Vec2& pos) const{
        // For numerical stability we add -eps to the upper bound of the clamping, to avoid the case when we are exactly on the edge and try to access out of bounds data.
        float x = std::clamp(pos.x, 0.0f, static_cast<float>(n_cols - 2) - eps);
        float y = std::clamp(pos.y, 0.0f, static_cast<float>(n_rows - 2) - eps);

        size_t col = static_cast<size_t>(x);
        size_t row = static_cast<size_t>(y);

        float x_frac = x - static_cast<float>(col);
        float y_frac = y - static_cast<float>(row);

        float v00 = at(row, col);
        float v01 = at(row, col + 1);
        float v10 = at(row + 1, col);
        float v11 = at(row + 1, col + 1);

        return bi_interpolate(v00, v01, v10, v11, x_frac, y_frac);
    }


    void min_max_normalize() 
    {
        if (value.empty()) return;

        float min_val = *std::min_element(value.begin(), value.end());
        float max_val = *std::max_element(value.begin(), value.end());

        if (max_val == min_val) {
            std::fill(value.begin(), value.end(), 0.0f);
            return;
        }

        for (auto& val : value) {
            val = (val - min_val) * 255.0f / (max_val - min_val);
        }
    }
};

struct Particle{
    Vec2 position;
    Vec2 velocity;
};

struct HashPoint{
    Vec2 pos;
    uint lineID;
};

class SpatialHash{
    private:
        std::unordered_map<long long, std::vector<HashPoint>> grid;
        float d_sep;    // Separation of the grid cells
        float d_test;   // Test distance for checking proximity = 0.5 d_sep

    public:
        SpatialHash(float separation) : d_sep(separation*2.0f) {d_test = 0.9*separation;}

        long long getKey(const Vec2& pos) const {
            int x_index = static_cast<int>(std::floor(pos.x / d_sep));
            int y_index = static_cast<int>(std::floor(pos.y / d_sep));
            return (static_cast<long long>(x_index) << 32) | static_cast<unsigned int>(y_index);
        }

        long long makeKey(int ix, int iy) const {
            return (static_cast<long long>(ix) << 32) | static_cast<unsigned int>(iy);
        }

        void insert(const Vec2& pos, uint lineID) {
            long long key = getKey(pos);
            grid[key].push_back({pos, lineID});
        }

        bool isTooClose(const Vec2& pos, uint currentLineID) const 
        {
            int ix = static_cast<int>(std::floor(pos.x / d_sep));
            int iy = static_cast<int>(std::floor(pos.y / d_sep));

            // Check the 3x3 block of cells around the point
            for (int i = ix - 1; i <= ix + 1; ++i) 
            {
                for (int j = iy - 1; j <= iy + 1; ++j) 
                {
                    long long key = makeKey(i, j);

                    auto it = grid.find(key);
                    if (it != grid.end()) 
                    {
                        //std::cout << "Checking cell (" << i << ", " << j << ") with " << it->second.size() << " points." << std::endl;
                        for (const auto& point : it->second) 
                        {
                            
                            // We ignore current line
                            if (point.lineID == currentLineID)
                                continue;

                            // Calculate squared distance -> sqrt has lower performance
                            float dx = point.pos.x - pos.x;
                            float dy = point.pos.y - pos.y;
                            if ((dx * dx + dy * dy) < (d_test * d_test)) {
                                return true; 
                            }
                        }
                    }
                }
            }
            
            // No points are too close
            return false;
        }   
};