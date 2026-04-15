#include "visualization.h"

#include <algorithm>

void createGSImage(const PixelPlane* plane, sf::Image* img) {
    img->resize({static_cast<unsigned>(plane->n_cols), static_cast<unsigned>(plane->n_rows)});

    for (size_t i = 0; i < plane->n_rows; i++) // Rows
    {
        for (size_t j = 0; j < plane->n_cols; j++) // Cols
        {
            // Use at(row, col)
            float val_float = plane->at(i, j);
            std::uint8_t val = static_cast<std::uint8_t>(std::clamp(val_float, 0.0f, 255.0f));
            
            // setPixel uses (x, y) which is (col, row)
            img->setPixel({static_cast<unsigned>(j), static_cast<unsigned>(i)}, sf::Color(val, val, val));
        }
    }
}