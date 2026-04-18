#include <iostream>
#include <vector>
#include <deque>
#include <unordered_map>
#include <string>
#include <random>
#include <cmath>
#include <algorithm>

// External Libraries
#include <SFML/Graphics.hpp>
#include <H5Cpp.h>

// Local includes
#include "vec2.h"
#include "data_structures.h"
#include "helpers.h"
#include "integrators.h"
#include "seeds.h"
#include "lic.h"
#include "visualization.h"

// Modify to your needs :)
const std::string isabel_file = "../data/isabel_2d.h5";
const std::string metsim_file = "../data/metsim1_2d.h5";

using namespace H5;

// Inferno colormap approximation (matplotlib inferno)
// Values from 0 to 255, RGB
std::vector<sf::Color> inferno_colormap = {
    sf::Color(0, 0, 4, 128),      // 0.0
    sf::Color(31, 12, 72, 128),   // 0.1
    sf::Color(85, 15, 109, 128),  // 0.2
    sf::Color(135, 23, 126, 128), // 0.3
    sf::Color(177, 42, 124, 128), // 0.4
    sf::Color(208, 70, 107, 128), // 0.5
    sf::Color(225, 104, 85, 128), // 0.6
    sf::Color(237, 142, 63, 128), // 0.7
    sf::Color(248, 184, 45, 128), // 0.8
    sf::Color(252, 231, 37, 128)  // 0.9
};

sf::Color getInfernoColor(float value) {
    // value should be normalized 0-1
    value = std::clamp(value, 0.0f, 1.0f);
    float index = value * 9.0f; // 10 colors, indices 0-9
    size_t idx = static_cast<size_t>(index);
    float frac = index - idx;
    if (idx >= 9) return inferno_colormap[9];
    sf::Color c1 = inferno_colormap[idx];
    sf::Color c2 = inferno_colormap[idx + 1];
    // Interpolate RGB
    uint8_t r = static_cast<uint8_t>(c1.r + frac * (c2.r - c1.r));
    uint8_t g = static_cast<uint8_t>(c1.g + frac * (c2.g - c1.g));
    uint8_t b = static_cast<uint8_t>(c1.b + frac * (c2.b - c1.b));
    return sf::Color(r, g, b, 200); // alpha 200 for lowered opacity
}

sf::ConvexShape createArrowHead(const Vec2& position, const Vec2& direction, float length, float width, const sf::Color& color) {
    Vec2 dir = direction.normalized();
    if (dir.length() < 1e-6f) {
        return sf::ConvexShape();
    }

    Vec2 perp{-dir.y, dir.x};
    Vec2 base = position - dir * length;

    sf::ConvexShape triangle;
    triangle.setPointCount(3);
    triangle.setPoint(0, sf::Vector2f(position.x, position.y));
    triangle.setPoint(1, sf::Vector2f(base.x + perp.x * width * 0.5f, base.y + perp.y * width * 0.5f));
    triangle.setPoint(2, sf::Vector2f(base.x - perp.x * width * 0.5f, base.y - perp.y * width * 0.5f));
    triangle.setFillColor(color);
    return triangle;
}

void createInfernoImage(const PixelPlane* plane, sf::Image* img) {
    img->resize({static_cast<unsigned>(plane->n_cols), static_cast<unsigned>(plane->n_rows)});

    // Find min and max for normalization
    float min_val = *std::min_element(plane->value.begin(), plane->value.end());
    float max_val = *std::max_element(plane->value.begin(), plane->value.end());
    float range = max_val - min_val;
    if (range == 0) range = 1.0f;

    for (size_t i = 0; i < plane->n_rows; i++) // Rows
    {
        for (size_t j = 0; j < plane->n_cols; j++) // Cols
        {
            float val = plane->at(i, j);
            float normalized = (val - min_val) / range;
            sf::Color color = getInfernoColor(normalized);
            // setPixel uses (x, y) which is (col, row)
            img->setPixel({static_cast<unsigned>(j), static_cast<unsigned>(i)}, color);
        }
    }
}

int main(){
    std::cout << "Starting even-spaced field lines with inferno heatmap..." << std::endl;

    RawData data;

    // Read file
    readH5File(isabel_file, &data, true);

    // Generate evenly spaced field lines
    //std::vector<Line> field_lines = getEvenSeed(&data, 2.0f, 0.2f, 600);
    std::vector<Line> field_lines = getEvenSeed(&data, 6.0f, 0.5f, 200);

    // Create PixelPlane for velocity magnitude
    PixelPlane magnitude_plane;
    magnitude_plane.n_rows = data.n_rows;
    magnitude_plane.n_cols = data.n_cols;
    magnitude_plane.value.resize(data.n_rows * data.n_cols);

    for (size_t i = 0; i < data.n_rows; ++i) {
        for (size_t j = 0; j < data.n_cols; ++j) {
            Vec2 vel = data.at(i, j);
            float mag = std::sqrt(vel.x * vel.x + vel.y * vel.y);
            magnitude_plane.set(i, j, mag);
        }
    }

    // Create inferno image
    sf::Image heatmap_image;
    createInfernoImage(&magnitude_plane, &heatmap_image);

    // Create texture and sprite for heatmap
    sf::Texture heatmap_texture;
    if (!heatmap_texture.loadFromImage(heatmap_image)) {
        std::cerr << "Failed to load heatmap texture from image" << std::endl;
        return EXIT_FAILURE;
    }
    sf::Sprite heatmap_sprite(heatmap_texture);

    // Test of SFML window
    // ===================================================================================================

    // create the window
    const size_t window_width = 1000;
    const size_t window_height = 1000;
    size_t cols = data.n_cols;
    size_t rows = data.n_rows;

    sf::RenderWindow window(sf::VideoMode({window_width, window_height}), "Even Field Lines with Inferno Heatmap");

    // run the program as long as the window is open
    while (window.isOpen())
    {
        // check all the window's events that were triggered since the last iteration of the loop
        while (const std::optional event = window.pollEvent())
        {
            // "close requested" event: we close the window
            if (event->is<sf::Event::Closed>())
                window.close();
        }

        // clear the window with white color
        window.clear(sf::Color::White);

        // Setting up margins to avoid drawing on the edge
        float margin = 20.0f;  // pixels
        float world_margin = margin * static_cast<float>(cols - 1) / window_width;

        // Let the library handle mapping from world to screen coordinates
        sf::View gridView;
        gridView.setSize({cols - 1.0f - 2*world_margin, rows - 1.0f - 2*world_margin});
        gridView.setCenter({cols / 2.0f, rows / 2.0f});
        window.setView(gridView);

        // Draw the heatmap
        window.draw(heatmap_sprite);

        // Randomize if the line will have arrow
        size_t count = 0;

        // We process all lines and particles inside of them and draw them as line strips
        for (const auto& line : field_lines)
        {
            sf::VertexArray strip(sf::PrimitiveType::LineStrip);
            for (const auto& particle : line.line)
            {
                strip.append(sf::Vertex{{particle.position.x, particle.position.y}, sf::Color(25, 25, 25, 255)});
            }
            window.draw(strip);

            if (!line.line.empty()) {
                const Particle& tail = line.line.back();
                Vec2 direction = data.interpolate(tail.position);
                if (direction.length() < 1e-6f && line.line.size() >= 2) {
                    const Particle& prev = line.line[line.line.size() - 2];
                    direction = (tail.position - prev.position).normalized();
                }

               
                if (direction.length() >= 1e-6f && count%5 == 0) {
                    float arrow_length = 1.8f*1.7;
                    float arrow_width = 1.2f*1.7;
                    sf::ConvexShape arrow = createArrowHead(tail.position, direction, arrow_length, arrow_width, sf::Color(25, 25, 25, 255));
                    if (arrow.getPointCount() == 3) {
                        window.draw(arrow);
                    }
                }

                count++;
            }
        }

        // Display everything on the window
        window.display();
    }

    return 0;
}