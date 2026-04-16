#include <iostream>
#include <vector>
#include <deque>
#include <unordered_map>
#include <string>
#include <random>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <sstream>

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

struct RandomParams {
    size_t max_steps;
    float particle_percentage;
};

struct UniformParams {
    size_t max_steps;
    size_t cell_distance;
};

struct EvenParams {
    size_t max_steps;
    float line_distance;
};

int main() {
    std::cout << "Starting field lines method comparison..." << std::endl;

    RawData data;
    std::string name = "metsim";
    readH5File(metsim_file, &data);

    const float step_size = 0.1f; // constant

    // Different parameter sets for each method
    /* 
    // Isabel dataset settings
    const std::vector<RandomParams> random_params = {
        {100, 0.001f},
        {200, 0.001f},
        {500, 0.001f},
        {800, 0.001f}
    };
    const std::vector<UniformParams> uniform_params = {
        {100, 25},
        {200, 25},
        {500, 25},
        {800, 25}
    };
    const std::vector<EvenParams> even_params = {
        {100, 5.0f},
        {200, 5.0f},
        {500, 5.0f},
        {800, 5.0f}
    };
    */

    // Metsim dataset settings
    const std::vector<RandomParams> random_params = {
        {50, 0.01f},
        {100, 0.01f},
        {200, 0.01f},
        {500, 0.01f}
    };
    const std::vector<UniformParams> uniform_params = {
        {50, 10},
        {100, 10},
        {200, 10},
        {500, 10}
    };
    const std::vector<EvenParams> even_params = {
        {50, 3.0f},
        {100, 3.0f},
        {200, 3.0f},
        {500, 3.0f}
    };

    const std::filesystem::path output_dir = "../exercises/field_lines/img";
    std::filesystem::create_directories(output_dir);

    const size_t window_width = 1000;
    const size_t window_height = 1000;
    size_t cols = data.n_cols;
    size_t rows = data.n_rows;

    const float margin = 20.0f;
    const float world_margin = margin * static_cast<float>(cols - 1) / window_width;

    sf::View gridView;
    gridView.setSize({cols - 1.0f - 2*world_margin, rows - 1.0f - 2*world_margin});
    gridView.setCenter({cols / 2.0f, rows / 2.0f});

    auto renderAndSave = [&](const std::vector<Line>& field_lines, const std::string& filename) {
        try {
            sf::RenderTexture renderTexture(sf::Vector2u{static_cast<unsigned int>(window_width), static_cast<unsigned int>(window_height)});
            renderTexture.clear(sf::Color::White);
            renderTexture.setView(gridView);

            for (const auto& line : field_lines) {
                sf::VertexArray strip(sf::PrimitiveType::LineStrip);
                for (const auto& particle : line.line) {
                    strip.append(sf::Vertex{{particle.position.x, particle.position.y}, sf::Color::Black});
                }
                renderTexture.draw(strip);
            }

            renderTexture.display();
            const sf::Texture& texture = renderTexture.getTexture();
            sf::Image image = texture.copyToImage();
            return image.saveToFile((output_dir / filename).string());
        } catch (const sf::Exception& e) {
            std::cerr << "SFML exception: " << e.what() << std::endl;
            return false;
        }
    };

    // Method 1: Random seed + RK4
    for (const auto& params : random_params) {
        std::cout << "Random + RK4: max_steps=" << params.max_steps << ", particle_percentage=" << params.particle_percentage << std::endl;
        std::vector<Particle> seeds = getRandomSeed(&data, params.particle_percentage);
        std::vector<Line> field_lines = rk4Integrator(&seeds, &data, step_size, params.max_steps);
        std::string filename = name + "_random_rk4_max" + std::to_string(params.max_steps) + "_pct" + std::to_string(static_cast<int>(params.particle_percentage * 100)) + ".png";
        if (renderAndSave(field_lines, filename)) {
            std::cout << "  Saved " << filename << std::endl;
        }
    }

    // Method 2: Uniform seed + RK4
    for (const auto& params : uniform_params) {
        std::cout << "Uniform + RK4: max_steps=" << params.max_steps << ", cell_distance=" << params.cell_distance << std::endl;
        std::vector<Particle> seeds = getUniformSeed(&data, params.cell_distance);
        std::vector<Line> field_lines = rk4Integrator(&seeds, &data, step_size, params.max_steps);
        std::string filename = name + "_uniform_rk4_max" + std::to_string(params.max_steps) + "_cell" + std::to_string(params.cell_distance) + ".png";
        if (renderAndSave(field_lines, filename)) {
            std::cout << "  Saved " << filename << std::endl;
        }
    }

    // Method 3: Even seed
    for (const auto& params : even_params) {
        std::cout << "Even seed: max_steps=" << params.max_steps << ", line_distance=" << params.line_distance << std::endl;
        std::vector<Line> field_lines = getEvenSeed(&data, params.line_distance, step_size, params.max_steps);
        std::string filename = name + "_even_seed_max" + std::to_string(params.max_steps) + "_line" + std::to_string(static_cast<int>(params.line_distance)) + ".png";
        if (renderAndSave(field_lines, filename)) {
            std::cout << "  Saved " << filename << std::endl;
        }
    }

    return 0;
}
