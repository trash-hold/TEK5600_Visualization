#include <iostream>
#include <vector>
#include <string>
#include <random> 
#include <cmath>
#include <H5Cpp.h>

// Modify to your needs :)
const std::string isabel_file = "./data/isabel_2d.h5";
const std::string metsim_file = "./data/metsim_2d.h5";

using namespace H5;

// A struct to hold the raw data read from the files
struct RawData{
    size_t n_rows, n_cols; // dimensions of the data;
    std::vector<float> x;
    std::vector<float> y;
};

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
    Vec2 operator*(float scalar) const {
        return Vec2(x * scalar, y * scalar);
    }
    Vec2 operator+=(const Vec2& other) {
        x += other.x;
        y += other.y;
        return *this;
    }

    // Helper functions
    float length() const {
        return std::sqrt(x * x + y * y);
    }

    Vec2 normalized() {
        float len = length();
         return len > 0 ? *this * (1.0f / len) : *this;
    }
};

struct Particle{
    Vec2 position;
    Vec2 velocity;
};


struct GridField{
    size_t nx, ny; // dimensions of the grid
    std::vector<float> x; // x component of the field
    std::vector<float> y; // y component of the field

    // Member functions
    Vec2 at(size_t i, size_t j) const {
        size_t index = i*nx + j;
        return Vec2(x[index], y[index]);
    }
};

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
        size_t col = static_cast<size_t>(rand_col);
        size_t row = static_cast<size_t>(rand_row);

        // Calculate indices for the four surrounding points in the data
        size_t i_00 = row * data->n_cols + col;
        size_t i_01 = row * data->n_cols + (col + 1);
        size_t i_10 = (row + 1) * data->n_cols + col;
        size_t i_11 = (row + 1) * data->n_cols + (col + 1);

        float x_frac = rand_col - col;
        float y_frac = rand_row - row;

        float inter_x = bi_interpolate(data->x[i_00], data->x[i_01], data->x[i_10], data->x[i_11], x_frac, y_frac);
        float inter_y = bi_interpolate(data->y[i_00], data->y[i_01], data->y[i_10], data->y[i_11], x_frac, y_frac);

        // Assign the position and velocity to the seed particle
        Particle p;
        p.position = Vec2(rand_col, rand_row);
        p.velocity = Vec2(inter_x, inter_y);
        seed.push_back(p);
    }

    return seed;
}

std::vector<Particle> getUniformSeed(RawData* const data, size_t cell_distance)
{
    // Initialize the particle vector
    std::vector<Particle> seed;
    size_t est_size = ((data->n_rows / cell_distance) + 1) * ((data->n_cols / cell_distance) + 1)
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
        for(size_t col=0; col < n_cols; j+=cell_distance)
        {
            rand_col = dis(gen);
            rand_row = dis(gen);
            
            // Calculate indices for the four surrounding points in the data
            size_t i_00 = row*(data->n_cols) + col;
            size_t i_01 = row*(data->n_cols) + (col+1);
            size_t i_10 = (row+1)*(data->n_cols) + col;
            size_t i_11 = (row+1)*(data->n_cols) + (col+1);

            float inter_x = bi_interpolate(data->x[i_00], data->x[i_01], data->x[i_10], data->x[i_11], rand_col, rand_row);

            float inter_y = bi_interpolate(data->y[i_00], data->y[i_01], data->y[i_10], data->y[i_11], rand_col, rand_row);

            // Assign the position and velocity to the seed particle
            Particle p;
            p.position = Vec2(col + rand_col, row + rand_row);
            p.velocity = Vec2(inter_x, inter_y);
            seed.push_back(p);
        }
    }

    return seed;
}

int main(){
    std::cout << "Starting processing of the files..." << std::endl;

    RawData isabel_data;
    RawData metsim_data;

    readH5File(isabel_file, &isabel_data);


    return 0;
}