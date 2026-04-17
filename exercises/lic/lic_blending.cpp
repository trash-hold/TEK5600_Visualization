#include <iostream>
#include <vector>
#include <string>
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

// Blending mode enumeration
enum class BlendMode {
    Multiply,
    Screen,
    SoftLight
};

// Inferno colormap approximation (matplotlib inferno)
std::vector<sf::Color> inferno_colormap = {
    sf::Color(0, 0, 4, 255),      // 0.0
    sf::Color(31, 12, 72, 255),   // 0.1
    sf::Color(85, 15, 109, 255),  // 0.2
    sf::Color(135, 23, 126, 255), // 0.3
    sf::Color(177, 42, 124, 255), // 0.4
    sf::Color(208, 70, 107, 255), // 0.5
    sf::Color(225, 104, 85, 255), // 0.6
    sf::Color(237, 142, 63, 255), // 0.7
    sf::Color(248, 184, 45, 255), // 0.8
    sf::Color(252, 231, 37, 255)  // 0.9
};

sf::Color getInfernoColor(float value) {
    value = std::clamp(value, 0.0f, 1.0f);
    float index = value * 9.0f;
    size_t idx = static_cast<size_t>(index);
    float frac = index - idx;
    if (idx >= 9) return inferno_colormap[9];
    sf::Color c1 = inferno_colormap[idx];
    sf::Color c2 = inferno_colormap[idx + 1];
    uint8_t r = static_cast<uint8_t>(c1.r + frac * (c2.r - c1.r));
    uint8_t g = static_cast<uint8_t>(c1.g + frac * (c2.g - c1.g));
    uint8_t b = static_cast<uint8_t>(c1.b + frac * (c2.b - c1.b));
    return sf::Color(r, g, b, 255);
}

void createInfernoImage(const PixelPlane* plane, sf::Image* img) {
    img->resize({static_cast<unsigned>(plane->n_cols), static_cast<unsigned>(plane->n_rows)});

    float min_val = *std::min_element(plane->value.begin(), plane->value.end());
    float max_val = *std::max_element(plane->value.begin(), plane->value.end());
    float range = max_val - min_val;
    if (range == 0) range = 1.0f;

    for (size_t i = 0; i < plane->n_rows; i++) {
        for (size_t j = 0; j < plane->n_cols; j++) {
            float val = plane->at(i, j);
            float normalized = (val - min_val) / range;
            sf::Color color = getInfernoColor(normalized);
            img->setPixel({static_cast<unsigned>(j), static_cast<unsigned>(i)}, color);
        }
    }
}

// Blending functions - normalize colors to [0, 1]
sf::Color blendMultiply(const sf::Color& licPixel, const sf::Color& heatmapPixel) {
    float lr = licPixel.r / 255.0f;
    float lg = licPixel.g / 255.0f;
    float lb = licPixel.b / 255.0f;
    
    float hr = heatmapPixel.r / 255.0f;
    float hg = heatmapPixel.g / 255.0f;
    float hb = heatmapPixel.b / 255.0f;
    
    uint8_t r = static_cast<uint8_t>(lr * hr * 255.0f);
    uint8_t g = static_cast<uint8_t>(lg * hg * 255.0f);
    uint8_t b = static_cast<uint8_t>(lb * hb * 255.0f);
    
    return sf::Color(r, g, b, 255);
}

sf::Color blendScreen(const sf::Color& licPixel, const sf::Color& heatmapPixel) {
    float lr = licPixel.r / 255.0f;
    float lg = licPixel.g / 255.0f;
    float lb = licPixel.b / 255.0f;
    
    float hr = heatmapPixel.r / 255.0f;
    float hg = heatmapPixel.g / 255.0f;
    float hb = heatmapPixel.b / 255.0f;
    
    uint8_t r = static_cast<uint8_t>((1.0f - (1.0f - lr) * (1.0f - hr)) * 255.0f);
    uint8_t g = static_cast<uint8_t>((1.0f - (1.0f - lg) * (1.0f - hg)) * 255.0f);
    uint8_t b = static_cast<uint8_t>((1.0f - (1.0f - lb) * (1.0f - hb)) * 255.0f);
    
    return sf::Color(r, g, b, 255);
}

sf::Color blendSoftLight(const sf::Color& licPixel, const sf::Color& heatmapPixel) {
    float lr = licPixel.r / 255.0f;
    float lg = licPixel.g / 255.0f;
    float lb = licPixel.b / 255.0f;
    
    float hr = heatmapPixel.r / 255.0f;
    float hg = heatmapPixel.g / 255.0f;
    float hb = heatmapPixel.b / 255.0f;
    
    auto softLightComponent = [](float base, float blend) -> float {
        if (blend < 0.5f) {
            return base - (1.0f - 2.0f * blend) * base * (1.0f - base);
        } else {
            float g = (base < 0.25f) ? ((16.0f * base - 12.0f) * base + 4.0f) * base 
                                      : std::sqrt(base);
            return base + (2.0f * blend - 1.0f) * (g - base);
        }
    };
    
    uint8_t r = static_cast<uint8_t>(std::clamp(softLightComponent(lr, hr), 0.0f, 1.0f) * 255.0f);
    uint8_t g = static_cast<uint8_t>(std::clamp(softLightComponent(lg, hg), 0.0f, 1.0f) * 255.0f);
    uint8_t b = static_cast<uint8_t>(std::clamp(softLightComponent(lb, hb), 0.0f, 1.0f) * 255.0f);
    
    return sf::Color(r, g, b, 255);
}

sf::Image blendImages(const sf::Image& licImage, const sf::Image& heatmapImage, BlendMode mode) {
    sf::Image result;
    result.resize({licImage.getSize().x, licImage.getSize().y});
    
    for (unsigned i = 0; i < licImage.getSize().y; ++i) {
        for (unsigned j = 0; j < licImage.getSize().x; ++j) {
            sf::Color licPixel = licImage.getPixel({j, i});
            sf::Color heatmapPixel = heatmapImage.getPixel({j, i});
            
            sf::Color blended;
            switch (mode) {
                case BlendMode::Multiply:
                    blended = blendMultiply(licPixel, heatmapPixel);
                    break;
                case BlendMode::Screen:
                    blended = blendScreen(licPixel, heatmapPixel);
                    break;
                case BlendMode::SoftLight:
                    blended = blendSoftLight(licPixel, heatmapPixel);
                    break;
            }
            result.setPixel({j, i}, blended);
        }
    }
    
    return result;
}

sf::Image applyVelocityMask(const sf::Image& blendedImage, const RawData* data, size_t scale_factor) {
    sf::Image masked;
    masked.resize({blendedImage.getSize().x, blendedImage.getSize().y});
    
    constexpr float zero_threshold = 1e-4f;
    
    for (unsigned i = 0; i < blendedImage.getSize().y; ++i) {
        for (unsigned j = 0; j < blendedImage.getSize().x; ++j) {
            // Map scaled coordinates back to original data space
            Vec2 pos((j + 0.5f) / scale_factor, (i + 0.5f) / scale_factor);
            Vec2 vel = data->interpolate(pos);
            float mag = std::sqrt(vel.x * vel.x + vel.y * vel.y);
            
            sf::Color blended = blendedImage.getPixel({j, i});
            if (mag < zero_threshold) {
                masked.setPixel({j, i}, sf::Color(0, 0, 0, 255));
            } else {
                masked.setPixel({j, i}, blended);
            }
        }
    }
    
    return masked;
}

std::string getModeString(BlendMode mode) {
    switch (mode) {
        case BlendMode::Multiply: return "Multiply";
        case BlendMode::Screen: return "Screen";
        case BlendMode::SoftLight: return "Soft Light";
        default: return "Unknown";
    }
}

int main() {
    std::cout << "Starting LIC with Blending Modes..." << std::endl;

    RawData data;
    readH5File(metsim_file, &data);

    std::cout << data.n_rows << " x " << data.n_cols << std::endl;

    // Compute LIC
    std::cout << "Computing LIC..." << std::endl;
    
    
    //For Metasim
    PixelPlane lic_result = simpleScaledLIC(&data, 5, 0.1f, 30, rk4Integrator);

    // For Isabel
    //PixelPlane lic_result = simpleScaledLIC(&data, 2, 0.5f, 10, rk4Integrator);

    // Create velocity magnitude heatmap (scaled to match LIC dimensions)
    std::cout << "Creating heatmap..." << std::endl;
    const size_t scale_factor = 5;
    PixelPlane magnitude_plane;
    magnitude_plane.n_rows = lic_result.n_rows;
    magnitude_plane.n_cols = lic_result.n_cols;
    magnitude_plane.value.resize(magnitude_plane.n_rows * magnitude_plane.n_cols);

    for (size_t i = 0; i < magnitude_plane.n_rows; ++i) {
        for (size_t j = 0; j < magnitude_plane.n_cols; ++j) {
            // Map scaled coordinates back to original data coordinates
            Vec2 pos((j + 0.5f) / scale_factor, (i + 0.5f) / scale_factor);
            Vec2 vel = data.interpolate(pos);
            float mag = std::sqrt(vel.x * vel.x + vel.y * vel.y);
            magnitude_plane.set(i, j, mag);
        }
    }

    // Create images
    sf::Image lic_image;
    createGSImage(&lic_result, &lic_image);

    sf::Image heatmap_image;
    createInfernoImage(&magnitude_plane, &heatmap_image);

    // Textures and sprites
    sf::Texture lic_texture;
    if (!lic_texture.loadFromImage(lic_image)) {
        std::cerr << "Failed to load LIC texture\n";
        return 1;
    }
    sf::Sprite lic_sprite(lic_texture);

    // Window setup
    const size_t window_width = 1000;
    const size_t window_height = 1000;
    size_t cols = lic_result.n_cols;
    size_t rows = lic_result.n_rows;

    sf::RenderWindow window(sf::VideoMode({window_width, window_height}), "LIC with Blending Modes");

    // Current blend mode
    BlendMode current_mode = BlendMode::SoftLight;

    // Main loop
    while (window.isOpen()) {
        while (const std::optional event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>()) {
                window.close();
            } else if (const auto* key_event = event->getIf<sf::Event::KeyPressed>()) {
                if (key_event->code == sf::Keyboard::Key::Num1) {
                    current_mode = BlendMode::Multiply;
                    std::cout << "Switched to: Multiply" << std::endl;
                } else if (key_event->code == sf::Keyboard::Key::Num2) {
                    current_mode = BlendMode::Screen;
                    std::cout << "Switched to: Screen" << std::endl;
                } else if (key_event->code == sf::Keyboard::Key::Num3) {
                    current_mode = BlendMode::SoftLight;
                    std::cout << "Switched to: Soft Light" << std::endl;
                }
            }
        }

        window.clear(sf::Color::White);

        // View setup
        float margin = 20.0f;
        float world_margin = margin * static_cast<float>(cols - 1) / window_width;

        sf::View gridView;
        gridView.setSize({static_cast<float>(cols - 1) - 2 * world_margin, static_cast<float>(rows - 1) - 2 * world_margin});
        gridView.setCenter({cols / 2.0f, rows / 2.0f});
        window.setView(gridView);

        // Blend and display
        sf::Image blended = blendImages(lic_image, heatmap_image, current_mode);
        
        // Apply velocity mask to hide zero-velocity regions
        sf::Image masked = applyVelocityMask(blended, &data, scale_factor);
        
        sf::Texture blended_texture;
        if (!blended_texture.loadFromImage(masked)) {
            std::cerr << "Failed to load blended texture\n";
        }
        sf::Sprite blended_sprite(blended_texture);

        window.draw(blended_sprite);

        // Display mode in window title
        window.setTitle("LIC with Blending Modes - " + getModeString(current_mode) + " (1/2/3 to switch)");

        window.display();
    }

    return 0;
}
