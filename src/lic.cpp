#include "lic.h"

PixelPlane simpleLIC(const RawData* data, const float& step_size ,const size_t& kernel_length)
{
    // Preparing the noise plane with padding
    PixelPlane noise = {
        data->n_rows,
        data->n_cols,
        std::vector<float>(data->n_rows * data->n_cols) // We allocate the pixel plane with the padding included
    };

    // Fill the noise plane with random values
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0f, 1.0f); // Less quantized than 0-255 int 

    for (size_t i = 0; i < noise.value.size(); i++) {
        noise.value[i] = dis(gen);
    }

    // Preparing the output plane
    PixelPlane output = {
        data->n_rows,
        data->n_cols,
        std::vector<float>(data->n_rows * data->n_cols) // Output plane does not need padding, so we allocate it with the original dimensions
    };

    // For the simple LIC we do convolution for each pixel
    // By assuming the same size between RawData and output plane we get 1->1 mapping
    for (size_t i = 0; i < data->n_rows; ++i) {
        for (size_t j = 0; j < data->n_cols; ++j) {
            std::vector<Particle> seed = {Particle{Vec2(j+0.5f, i+0.5f), bi_interpolate(data, Vec2(j+0.5f, i+0.5f))}}; // We sample at the center of the pixel

            // The kernel length is interpreted as the radius of the convolution
            Line new_line = rk4Integrator(&seed, data, step_size, kernel_length)[0]; // We get the field line

            float sum = 0.0f;
            if (new_line.line.size() == 0) {
                output.set(i, j, noise.at(i, j)); // If the line is empty, we set the output to the noise value at that point
                continue;
            }

            for (auto& particle : new_line.line) 
            {
                sum += noise.interpolate(particle.position);
            }
            
            output.set(i, j, sum / new_line.line.size()); 
        }
    }

    output.min_max_normalize(); // We normalize the output to [0, 255] range for better visualization
    return output;
}