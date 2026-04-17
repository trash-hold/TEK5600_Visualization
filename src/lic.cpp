#include "lic.h"

template <typename F> 
PixelPlane simpleLIC(const RawData* data, const float& step_size, const size_t& kernel_length, F integrator)
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
    for (size_t i = 0; i < data->n_rows; i++) {
        for (size_t j = 0; j < data->n_cols; j++) {
            std::vector<Particle> seed = {Particle{Vec2(j+0.5f, i+0.5f), data->interpolate(Vec2(j+0.5f, i+0.5f))}}; // We sample at the center of the pixel

            // The kernel length is interpreted as the radius of the convolution
            Line new_line = integrator(&seed, data, step_size, kernel_length)[0]; // We get the field line

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

template PixelPlane simpleLIC(const RawData*, const float&, const size_t&, 
                             std::vector<Line> (*)(const std::vector<Particle> *seeds, const RawData *data, const float &step_size, const size_t &max_steps));

template <typename F2>
PixelPlane simpleScaledLIC(const RawData* data, const size_t &scale_factor, const float& step_size, const size_t& kernel_length, F2 integrator)
{
    if(scale_factor < 1)
        throw std::runtime_error("Scaling factor needs to be greater that 1");

    size_t rescaled_rows = data->n_rows * scale_factor;
    size_t rescaled_columns = data->n_cols * scale_factor;

    // Preparing the noise plane
    PixelPlane noise = {
        rescaled_rows,
        rescaled_columns,
        std::vector<float>(rescaled_rows * rescaled_columns),
        scale_factor
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
        rescaled_rows,
        rescaled_columns,
        std::vector<float>(rescaled_rows * rescaled_columns),
        scale_factor
    };

    // For the simple LIC we do convolution for each pixel
    for (size_t i = 0; i < rescaled_rows; i++) {
        for (size_t j = 0; j < rescaled_columns; j++) {
            Vec2 position = Vec2( (j+0.5)/scale_factor, (i+0.5)/scale_factor); 
            std::vector<Particle> seed = {Particle{position, data->interpolate(position)}}; // We sample at the center of the pixel

            // The kernel length is interpreted as the radius of the convolution
            Line new_line = integrator(&seed, data, step_size, kernel_length)[0]; // We get the field line

            float sum = 0.0f;
            if (new_line.line.size() == 0) {
                output.set(i, j, noise.at(i, j)); // If the line is empty, we set the output to the noise value at that point
                continue;
            }

            for (auto& particle : new_line.line) 
            {
                // We need to rescale
                sum += noise.interpolate(noise.scalePosition(particle.position));
            }
            
            output.set(i, j, sum / new_line.line.size()); 
        }
    }

    output.min_max_normalize(); // We normalize the output to [0, 255] range for better visualization
    return output;
}

template PixelPlane simpleScaledLIC(const RawData*, const size_t&, const float&, const size_t&, 
                             std::vector<Line> (*)(const std::vector<Particle> *seeds, const RawData *data, const float &step_size, const size_t &max_steps));