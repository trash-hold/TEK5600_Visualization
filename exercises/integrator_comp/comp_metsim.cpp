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

int main(){
    std::cout << "Starting field lines processing..." << std::endl;

    RawData data;

    // Read file
    readH5File(metsim_file, &data);

    // Generate seeds
    std::vector<Particle> seeds = getUniformSeed(&data, 5);

    std::cout << "Generated " << seeds.size() << " seeds" << std::endl;

    // Experiment configuration
    const std::vector<float> step_sizes = {0.1f, 0.2f, 0.5f, 1.0f, 2.0f};
    const float total_length = 20.0f;
    const std::filesystem::path output_dir = "../exercises/integrator_comp/img";
    std::filesystem::create_directories(output_dir);

    const size_t window_width = 1000;
    const size_t window_height = 1000;
    size_t cols = data.n_cols;
    size_t rows = data.n_rows;

    const float margin = 20.0f;  // pixels
    const float world_margin = margin * static_cast<float>(cols - 1) / window_width;

    sf::View gridView;
    gridView.setSize({cols - 1.0f - 2*world_margin, rows - 1.0f - 2*world_margin});
    gridView.setCenter({cols / 2.0f, rows / 2.0f});

    auto formatStepSize = [&](float value) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(2) << value;
        return ss.str();
    };

    auto renderAndSave = [&](const std::vector<Line>& field_lines, const std::string& filename) {
        try {
            sf::RenderTexture renderTexture(sf::Vector2u{static_cast<unsigned int>(window_width), static_cast<unsigned int>(window_height)});
            renderTexture.clear(sf::Color::White);
            renderTexture.setView(gridView);

            for (const auto& line : field_lines)
            {
                sf::VertexArray strip(sf::PrimitiveType::LineStrip);
                for (const auto& particle : line.line)
                {
                    strip.append(sf::Vertex{{particle.position.x, particle.position.y}, sf::Color::Black});
                }
                renderTexture.draw(strip);
            }

            renderTexture.display();
            const sf::Texture& texture = renderTexture.getTexture();
            sf::Image image = texture.copyToImage();
            return image.saveToFile((output_dir / filename).string());
        }
        catch (const sf::Exception& e) {
            std::cerr << "SFML exception: " << e.what() << std::endl;
            return false;
        }
    };

    for (float step_size : step_sizes)
    {
        const float exact_steps = total_length / step_size;
        const size_t max_steps = static_cast<size_t>(std::round(exact_steps));

        if (std::fabs(exact_steps - static_cast<float>(max_steps)) > 1e-6f)
        {
            std::cerr << "step_size " << step_size << " does not divide 10 exactly; skipping." << std::endl;
            continue;
        }

        std::cout << "Running experiment: step_size=" << step_size << " max_steps=" << max_steps << std::endl;

        struct IntegratorInfo { const char* name; std::vector<Line> (*func)(const std::vector<Particle>*, const RawData*, const float&, const size_t&); } integrators[] = {
            {"euler", eulerIntegrator},
            {"rk4", rk4Integrator}
        };

        for (const auto& integrator : integrators)
        {
            std::vector<Line> field_lines = integrator.func(&seeds, &data, step_size, max_steps);
            std::cout << "  " << integrator.name << " -> " << field_lines.size() << " lines" << std::endl;

            const std::string filename = "metsim_" + std::string(integrator.name) + "_step" + formatStepSize(step_size) + ".png";
            if (renderAndSave(field_lines, filename))
            {
                std::cout << "  Saved " << filename << std::endl;
            }
            else
            {
                std::cerr << "  Failed to save " << filename << std::endl;
            }
        }
    }

    return 0;
}