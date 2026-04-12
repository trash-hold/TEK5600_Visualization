#include <iostream>
#include <vector>
#include <deque>
#include <string>
#include <random> 
#include <cmath>

// External Libraries
#include <SFML/Graphics.hpp>
#include <H5Cpp.h>

// Modify to your needs :)
const std::string isabel_file = "../data/isabel_2d.h5";
const std::string metsim_file = "../data/metsim_2d.h5";

using namespace H5;

// ==================================================================================================
//  Primitives
// ==================================================================================================

// Primitive for the vector field
struct Vec2{
    float x, y;

    // Member Functions
    Vec2(float xx = 0.0f, float yy = 0.0f) : x(xx), y(yy) {}

    // Operator Overloads
    Vec2 operator+(const Vec2& other) const {
        return Vec2(x + other.x, y + other.y);
    }
    Vec2 operator-(const Vec2& other) const {
        return Vec2(x - other.x, y - other.y);
    }
    Vec2 operator-() const {
        return Vec2(-x, -y);
    }
    Vec2 operator*(float scalar) const {
        return Vec2(x * scalar, y * scalar);
    }
    Vec2 operator+=(const Vec2& other) {
        x += other.x;
        y += other.y;
        return *this;
    }
    Vec2 operator-=(const Vec2& other) {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    // Helper functions
    float length() const {
        return std::sqrt(x * x + y * y);
    }

    Vec2 normalized() const{
        float len = length();
        
        if (len > 0.0f) 
        {
            return Vec2(x / len, y / len);
        } 
        else 
        { 
            return Vec2(0.0f, 0.0f); // Return zero vector if length is zero to avoid division by zero
        }
    }
};


// A struct to hold the raw data read from the files
struct RawData{
    size_t n_rows, n_cols; // dimensions of the data;
    std::vector<float> x;
    std::vector<float> y;

    Vec2 at(size_t i, size_t j) const {
        if (i >= n_rows || j >= n_cols) {
            throw std::out_of_range("Index out of bounds");
        }
        size_t index = i*n_cols + j;
        return Vec2(x[index], y[index]);
    }
};


struct Particle{
    Vec2 position;
    Vec2 velocity;
};


// ==================================================================================================
//  Helpers
// ==================================================================================================

// Function to read the HDF5 file -- there is assumed structure of the file and equality of dimensions for both x and y components
void readH5File(const std::string& filename, RawData* data){
    try {
        // Open file in read-only mode
        H5File file(filename, H5F_ACC_RDONLY);

        // There are two datasets in both files:
        DataSet x_dataset = file.openDataSet("/Velocity/X-comp");
        DataSet y_dataset = file.openDataSet("/Velocity/Y-comp");

        // Getting size of the datasets (assuming both datasets have the same dimensions)
        hsize_t dims[2];
        DataSpace space = x_dataset.getSpace();
        space.getSimpleExtentDims(dims);

        // The read data will be stored as 1D vector -- it's faster
        data->n_rows = dims[0];
        data->n_cols = dims[1];
        data->x.resize(dims[0] * dims[1]);
        data->y.resize(dims[0] * dims[1]);

        x_dataset.read(data->x.data(), PredType::NATIVE_FLOAT);
        y_dataset.read(data->y.data(), PredType::NATIVE_FLOAT);

        std::cout << dims[0] << " x " << dims[1] << "\n";

        file.close();
    }
    catch (H5::Exception& e) {
        std::cerr << "Error reading file: " << e.getCDetailMsg() << std::endl;
    }
}

float bi_interpolate(float v00, float v01, float v10, float v11, float x_frac, float y_frac) {
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

Vec2 bi_interpolate(const RawData* data, const Vec2 &pos)
{
    // Calculate the indices of the four surrounding particles
    size_t cols = data->n_cols;
    size_t rows = data->n_rows;

    size_t col = static_cast<size_t>(pos.x);
    size_t row = static_cast<size_t>(pos.y);

    // Interpolate 
    float x_frac = pos.x - col;
    float y_frac = pos.y - row;

    // Data
    Vec2 v00 = data->at(row, col);
    Vec2 v01 = col + 1 < cols ? data->at(row, col + 1) : v00;
    Vec2 v10 = row + 1 < rows ? data->at(row + 1, col) : v00;
    Vec2 v11 = row + 1 < rows && col + 1 < cols ? data->at(row + 1, col + 1) : v00;

    float inter_x = bi_interpolate(v00.x, v01.x, v10.x, v11.x, x_frac, y_frac);
    float inter_y = bi_interpolate(v00.y, v01.y, v10.y, v11.y, x_frac, y_frac);

    return {Vec2(inter_x, inter_y)};
}


// ==================================================================================================
//  Seed strategies
// ==================================================================================================
std::vector<Particle> getRandomSeed(RawData* const data, float particle_percentage)
{
    size_t num_particles = static_cast<size_t>(data->n_rows * data->n_cols * particle_percentage);

    // Initialize the particle vector
    std::vector<Particle> seed;
    seed.reserve(num_particles); // reserve space for all particles
    
    // Initialize random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_x(0.0f, static_cast<float>(data->n_cols - 1)); 
    std::uniform_real_distribution<> dis_y(0.0f, static_cast<float>(data->n_rows - 1));

    for (size_t i = 0; i < num_particles; ++i)
    {
        float rand_col = dis_x(gen);
        float rand_row = dis_y(gen);

        // Assign the position and velocity to the seed particle
        Particle p;
        p.position = Vec2(rand_col, rand_row);
        p.velocity = bi_interpolate(data, Vec2(rand_col, rand_row));
        seed.push_back(p);
    }

    return seed;
}

std::vector<Particle> getUniformSeed(RawData* const data, size_t cell_distance)
{
    // Initialize the particle vector
    std::vector<Particle> seed;
    size_t est_size = ((data->n_rows / cell_distance) + 1) * ((data->n_cols / cell_distance) + 1);
    seed.reserve(est_size); // reserve space for all particles

    // Initialize random number generator -- to nuance the positioning of the seed particles within the grid cells
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0f, 1.0f);

    // Allocate before loop
    float rand_col, rand_row;
    size_t n_rows = data->n_rows - 1; // To avoid out of bounds access when interpolating
    size_t n_cols = data->n_cols - 1;

    // Now we pick seed points and assign values based on interpolation
    for (size_t row=0; row < n_rows; row+=cell_distance)
    {
        for(size_t col=0; col < n_cols; col+=cell_distance)
        {
            rand_col = dis(gen);
            rand_row = dis(gen);
        
            // Getting surrounding points for interpolation
            Vec2 v00 = data->at(row, col);
            Vec2 v01 = data->at(row, col + 1);
            Vec2 v10 = data->at(row + 1, col);
            Vec2 v11 = data->at(row + 1, col + 1);

            float inter_x = bi_interpolate(v00.x, v01.x, v10.x, v11.x, rand_col, rand_row);
            float inter_y = bi_interpolate(v00.y, v01.y, v10.y, v11.y, rand_col, rand_row);

            // Assign the position and velocity to the seed particle
            Particle p;
            p.position = Vec2(col + rand_col, row + rand_row);
            p.velocity = Vec2(inter_x, inter_y);
            seed.push_back(p);
        }
    }

    return seed;
}

// ==================================================================================================
//  Line Integrators
// ==================================================================================================

struct Line{
    std::vector<Particle> line;
};


std::vector<Line> eulerIntegrator(const std::vector<Particle>* seeds, const RawData* data, const float &step_size, const size_t &max_steps)
{
    std::vector<Line> field_lines; 
    field_lines.reserve(seeds->size()); // We compute field-lines for all generated seeds

    size_t n_rows = data->n_rows;
    size_t n_cols = data->n_cols;

    for (const auto& seed : *seeds)
    {
        // Prepare temporary deque for line object      
        std::deque<Particle> deq_line;
        
        // Forward pass
        Vec2 f_current_pos = seed.position;
        Vec2 f_current_v = seed.velocity;

        // Backward pass
        Vec2 b_current_pos = seed.position;
        Vec2 b_current_v = -seed.velocity;

        deq_line.push_back(seed);

        for(int i = 0; i < max_steps; i++)
        {
            // Check if the new position is out of bounds
            if (f_current_pos.x < 0 || f_current_pos.x >= n_cols || f_current_pos.y < 0 || f_current_pos.y >= n_rows) {
                break; // Stop integrating if we go out of bounds
            }

            if (b_current_pos.x < 0 || b_current_pos.x >= n_cols || b_current_pos.y < 0 || b_current_pos.y >= n_rows) {
                break; // Stop integrating if we go out of bounds
            }

            // Sample the velcocity at the current positions
            f_current_v = bi_interpolate(data, f_current_pos);
            b_current_v = bi_interpolate(data, b_current_pos);

            // Calculate the new positions using Euler's method
            f_current_pos += f_current_v.normalized() * step_size;
            b_current_pos -= b_current_v.normalized() * step_size;

            // Add the new points to the line
            deq_line.push_back(Particle{f_current_pos, f_current_v});
            deq_line.push_front(Particle{b_current_pos, b_current_v});

        }

        // Deque to vector
        std::vector<Particle> particles(deq_line.begin(), deq_line.end());
        field_lines.push_back(Line{std::move(particles)});
    }
    return field_lines;
}

// ==================================================================================================
//  Visualization
// ==================================================================================================



int main(){
    std::cout << "Starting processing of the files..." << std::endl;

    RawData isabel_data;
    RawData metsim_data;

    // Read file
    readH5File(isabel_file, &isabel_data);

    // Generate seeds
    std::vector<Particle> seeds = getUniformSeed(&isabel_data, 10);

    // Calculate field lines
    std::vector<Line> field_lines = eulerIntegrator(&seeds, &isabel_data, 0.1f, 100);

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