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

int main(){
    std::cout << "Starting processing of the files..." << std::endl;

    RawData isabel_data;
    RawData metsim_data;

    // Read file
    readH5File(metsim_file, &isabel_data);

    // Run LIC
    PixelPlane lic_result = simpleLIC(&isabel_data, 0.1f, 30);

    // Prepare objects for visualization
    sf::Image lic_image;
    createGSImage(&lic_result, &lic_image);

    // Easiest way to display is through texture + sprite
    sf::Texture lic_texture;

    if (!lic_texture.loadFromImage(lic_image)) {
        std::cerr << "Failed to load texture\n";
        return 1;
    }

    sf::Sprite lic_sprite(lic_texture);

    // Test of SFML window
    // ===================================================================================================
    
    // create the window
    const size_t window_width = 1000;
    const size_t window_height = 1000;
    size_t cols = isabel_data.n_cols;
    size_t rows = isabel_data.n_rows;

    sf::RenderWindow window(sf::VideoMode({window_width, window_height}), "Test Window");

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
        
        // Setiting up margins to avoid drawing on the edge
        float margin = 20.0f;  // pixels
        float world_margin = margin * static_cast<float>(cols- 1) / window_width;

        // Let the library handle mapping from world to screen coordinates
        sf::View gridView;
        gridView.setSize({cols - 1.0f - 2*world_margin, rows - 1.0f - 2*world_margin});
        gridView.setCenter({cols / 2.0f, rows / 2.0f});
        window.setView(gridView); 

        // Display the spirte
        window.draw(lic_sprite);

        // Display everthing on the window
        window.display();
    }

    return 0;
}

/*
int main(){
    std::cout << "Starting processing of the files..." << std::endl;

    RawData isabel_data;
    RawData metsim_data;

    // Read file
    readH5File(isabel_file, &isabel_data);

    // Generate seeds
    std::vector<Particle> seeds = getUniformSeed(&isabel_data, 10);

    // Calculate field lines
    //std::vector<Line> field_lines = rk4Integrator(&seeds, &isabel_data, 0.1f, 100);
    std::vector<Line> field_lines = getEvenSeed(&isabel_data, 10.0f, 0.2f, 700);

    // Test of SFML window
    // ===================================================================================================
    
    // create the window
    const size_t window_width = 1000;
    const size_t window_height = 1000;
    size_t cols = isabel_data.n_cols;
    size_t rows = isabel_data.n_rows;

    sf::RenderWindow window(sf::VideoMode({window_width, window_height}), "Test Window");

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
        
        // Setiting up margins to avoid drawing on the edge
        float margin = 20.0f;  // pixels
        float world_margin = margin * static_cast<float>(cols- 1) / window_width;

        // Let the library handle mapping from world to screen coordinates
        sf::View gridView;
        gridView.setSize({cols - 1.0f - 2*world_margin, rows - 1.0f - 2*world_margin});
        gridView.setCenter({cols / 2.0f, rows / 2.0f});
        window.setView(gridView); 

        // We process all lines and particles inside of them and draw them as line strips
        for (const auto& line : field_lines)
        {   
            sf::VertexArray strip(sf::PrimitiveType::LineStrip);
            for (const auto& particle : line.line)
            {
                strip.append(sf::Vertex{{particle.position.x, particle.position.y}, sf::Color::Black});
            }
            window.draw(strip);
        }

        // Display everthing on the window
        window.display();
    }

    return 0;
}
*/