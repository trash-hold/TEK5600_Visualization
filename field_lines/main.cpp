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

// Modify to your needs :)
const std::string isabel_file = "../data/isabel_2d.h5";
const std::string metsim_file = "../data/metsim1_2d.h5";

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

Vec2 bi_interpolate(const RawData* data, const Vec2& pos)
{
    const float eps = 1e-4f;

    // For numerical stability we add -eps to the upper bound of the clamping, to avoid the case when we are exactly on the edge and try to access out of bounds data.
    float x = std::clamp(pos.x, 0.0f, static_cast<float>(data->n_cols - 2) - eps);
    float y = std::clamp(pos.y, 0.0f, static_cast<float>(data->n_rows - 2) - eps);

    size_t col = static_cast<size_t>(x);
    size_t row = static_cast<size_t>(y);

    float x_frac = x - static_cast<float>(col);
    float y_frac = y - static_cast<float>(row);

    Vec2 v00 = data->at(row, col);
    Vec2 v01 = data->at(row, col + 1);
    Vec2 v10 = data->at(row + 1, col);
    Vec2 v11 = data->at(row + 1, col + 1);

    float inter_x = bi_interpolate(v00.x, v01.x, v10.x, v11.x, x_frac, y_frac);
    float inter_y = bi_interpolate(v00.y, v01.y, v10.y, v11.y, x_frac, y_frac);

    return Vec2(inter_x, inter_y);
}

// ==================================================================================================
//  Line Integrators
// ==================================================================================================

struct Line{
    std::vector<Particle> line;
};


Particle euler_step(const Vec2& pos, const RawData* data, const float& step_size) {
    Vec2 velocity = bi_interpolate(data, pos).normalized() * step_size;
    Vec2 new_pos = pos + velocity;
    return Particle{new_pos, bi_interpolate(data, new_pos)};
}

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
        bool f_inbounds = true;

        // Backward pass
        Vec2 b_current_pos = seed.position;
        Vec2 b_current_v = -seed.velocity;
        bool b_inbounds = true;

        deq_line.push_back(seed);

        for(int i = 0; i < max_steps && (f_inbounds || b_inbounds); i++)
        {
            // Check if the new position is out of bounds
            if (f_current_pos.x < 0 || f_current_pos.x >= n_cols || f_current_pos.y < 0 || f_current_pos.y >= n_rows) {
                f_inbounds = false;
            }

            if (b_current_pos.x < 0 || b_current_pos.x >= n_cols || b_current_pos.y < 0 || b_current_pos.y >= n_rows) {
                b_inbounds = false;
            }
            
            if(f_inbounds)
            {
                Particle forward = euler_step(f_current_pos, data, step_size);
                f_current_pos = forward.position;
                // Add the new point to the line
                deq_line.push_back(forward);
            }

            if(b_inbounds)
            {
                Particle backward = euler_step(b_current_pos, data, -step_size);
                b_current_pos = backward.position;
                // Add the new point to the line
                deq_line.push_front(backward);
            }

        }

        // Deque to vector
        std::vector<Particle> particles(deq_line.begin(), deq_line.end());
        field_lines.push_back(Line{std::move(particles)});
    }
    return field_lines;
}

Particle rk4_step(const Vec2& pos, const RawData* data, const float& step_size) {
    Vec2 k1 = bi_interpolate(data, pos).normalized() * step_size;
    Vec2 k2 = bi_interpolate(data, pos + k1*0.5f).normalized() * step_size;
    Vec2 k3 = bi_interpolate(data, pos + k2*0.5f).normalized() * step_size;
    Vec2 k4 = bi_interpolate(data, pos + k3).normalized() * step_size;

    Vec2 new_pos = pos + (k1 + k2*2.0f + k3*2.0f + k4) * (1.0f / 6.0f);
    return Particle{new_pos, bi_interpolate(data, new_pos)};
}

std::vector<Line> rk4Integrator(const std::vector<Particle>* seeds, const RawData* data, const float &step_size, const size_t &max_steps)
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
        bool f_inbounds = true;

        // Backward pass
        Vec2 b_current_pos = seed.position;
        bool b_inbounds = true;

        deq_line.push_back(seed);

        for(int i = 0; i < max_steps && (f_inbounds || b_inbounds); i++)
        {
            // Check if the new position is out of bounds
            if (f_current_pos.x < 0 || f_current_pos.x >= n_cols || f_current_pos.y < 0 || f_current_pos.y >= n_rows) {
                f_inbounds = false;
            }

            if (b_current_pos.x < 0 || b_current_pos.x >= n_cols || b_current_pos.y < 0 || b_current_pos.y >= n_rows) {
                b_inbounds = false;
            }

        
            if(f_inbounds)
            {
                Particle forward = rk4_step(f_current_pos, data, step_size);
                f_current_pos = forward.position;
                deq_line.push_back(forward);
            }

            if(b_inbounds)
            {
                Particle backward = rk4_step(b_current_pos, data, -step_size);
                b_current_pos = backward.position;
                deq_line.push_front(backward);
            }

        }

        // Deque to vector
        std::vector<Particle> particles(deq_line.begin(), deq_line.end());
        field_lines.push_back(Line{std::move(particles)});
    }
    return field_lines;
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

bool growFieldLine(SpatialHash& hash, const RawData* data, std::vector<Line>& field_lines, const Vec2& seed, const uint& lineID, const float &step_size, const size_t &max_steps)
{
    // Check if seed is not too close
    if (hash.isTooClose(seed, lineID)) {
        return false; // Seed is too close to existing line, do not grow
    }

    std::deque<Particle> newLine;
    
    // Using RK4 integrator for better accuracy
    Vec2 f_current_pos = seed;
    Vec2 b_current_pos = seed;
    bool f_inbounds = true;
    bool b_inbounds = true; 
    bool ended_by_distance = false; // Flag to know that proimity violation was the reason for stopping the line growth, not the out of bounds. This is used to accept lines that are shorter than the nominal length, but we have no more space to grow them without violating the distance constraint.

    // Calculate the line in both directions until we hit the max steps or go out of bounds
    for(int i = 0; i < max_steps && (f_inbounds || b_inbounds); i++)
    {
        // Check if the new position is out of bounds
        // Forward step
        if (data->isValid(static_cast<size_t>(f_current_pos.y), static_cast<size_t>(f_current_pos.x))) 
        {
            // Standard step
            Particle forward = rk4_step(f_current_pos, data, step_size);
            f_current_pos = forward.position;
            
            // Checking if the new point is too close to existing lines
            if (hash.isTooClose(f_current_pos, lineID) == false) 
                newLine.push_back(forward);
            else
            {
                f_inbounds = false;
                ended_by_distance = true;
            }
                
            
        }
        else
        {
            f_inbounds = false;
        }

        // Backward step
        if (data->isValid(static_cast<size_t>(b_current_pos.y), static_cast<size_t>(b_current_pos.x))) 
        {
            Particle backward = rk4_step(b_current_pos, data, -step_size);
            b_current_pos = backward.position;

            if (hash.isTooClose(b_current_pos, lineID) == false) 
                newLine.push_front(backward);
            else
            {
                b_inbounds = false;
                ended_by_distance = true;
            }
        }
        else
        {
             b_inbounds = false;
        }

    }

    // The algorithm will have problems with longer lines, so we accept lines even if they don't have the "nominal" length
    size_t min_length = static_cast<size_t>(max_steps * 0.7f);
    if (min_length < 30) // We set a lower bound on the minimum length, to avoid accepting very short lines that are not informative
        min_length = 30;

    if (newLine.size() < min_length && ended_by_distance) 
    {
        return false; // Line is too short, do not add
    }

    // The line got computed without distance violations, now it's commited to the hash map
    std::vector<Particle> particles(newLine.begin(), newLine.end());
    for (const auto& particle : particles) 
    {
        //std::cout << "Inserting point (" << particle.position.x << ", " << particle.position.y << ") of line " << lineID << " into spatial hash." << std::endl;
        hash.insert(particle.position, lineID);
    }

    field_lines.push_back(Line{std::move(particles)});
    //std::cout << "Finished processing seeds. Total lines generated: " << field_lines.size() << std::endl;
    // Line successfully grown and added to the field lines
    return true;
}

std::vector<Line> getEvenSeed(RawData* const data, const float& line_distance, const float& step_size, const size_t& max_steps)
{
    // Initializing objects
    SpatialHash hash(line_distance);
    std::vector<Line> field_lines;
    uint lineID = 0;

    // Reserve space for field lines
    field_lines.reserve(static_cast<size_t>(data->n_rows * data->n_cols / line_distance));

    // Random seed generation for the first line
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_x(0.0f, static_cast<float>(data->n_cols - 1)); 
    std::uniform_real_distribution<> dis_y(0.0f, static_cast<float>(data->n_rows - 1));

    Particle start_seed;

    // We try max 1000 seed points
    for(int i=0; i < 1000; i++)
    {
        float rand_col = dis_x(gen);
        float rand_row = dis_y(gen);
        
        Vec2 velocity = bi_interpolate(data, Vec2(rand_col, rand_row));

        if (velocity.length() > 1e-3f) // We skip points with very low velocity to avoid drawing field lines in the areas where there is no flow
        {
            start_seed.position = Vec2(rand_col, rand_row);
            start_seed.velocity = velocity;

            break;
        }
    }

    // We check if the seed is valid
    if(start_seed.velocity.length() < 1e-3f)
    {
        std::cerr << "Failed to find a valid seed point after 1000 attempts." << std::endl;
        return field_lines; // Return empty vector
    }

    // Now we make a queue to keep track of unprocessed particles
    std::deque<Particle> seed_queue;
    seed_queue.push_back(start_seed);

    // To avoid seed explosion, only portion of the seeds are picked
    uint seed_skip = static_cast<uint>(line_distance * 0.05f / step_size);

    if (seed_skip < 1) {
        seed_skip = 1;
    }

    // We process the queue until it's empty
    while (!seed_queue.empty()) // We also set a lower bound on the number of lines to be generated, to avoid the case when the seed points are generated in the areas with no flow and we end up with no lines.
    {
        Particle current_seed = seed_queue.front();
        seed_queue.pop_front();

        // We try to grow a line from the seed
        if (growFieldLine(hash, data, field_lines, current_seed.position, lineID, step_size, max_steps))
        {
            // If the line was successfully grown, we add new seeds at a distance of line_distance from the line points
            const Line& new_line = field_lines.back();
            for (int i = 0; i < new_line.line.size(); i+=seed_skip)
            {

                Particle particle = new_line.line[i];
                Vec2 pos = particle.position;
                
                // We check if the velocity at the point is not too low, to avoid placing seeds in areas with no flow
                if (particle.velocity.length() < 1e-4f)
                    continue;

                // We add seeds in that are perpendicular to the velocity vector at the point
                Vec2 perp = Vec2(-particle.velocity.y, particle.velocity.x).normalized();
                std::vector<Vec2> candidate_seeds = {
                    pos + perp * line_distance, // Perpendicular in one direction
                    pos - perp * line_distance  // Perpendicular in the other direction
                };

                for (const auto& candidate : candidate_seeds) 
                {
                    if (data->isValid(static_cast<size_t>(candidate.y), static_cast<size_t>(candidate.x)) && hash.isTooClose(candidate, lineID) == false) 
                    {
                        seed_queue.push_back(Particle{candidate, bi_interpolate(data, candidate)});
                    }
                }
            }

            // Incrementing line ID for the next line
            lineID++;
        }
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