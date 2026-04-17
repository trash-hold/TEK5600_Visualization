#include <iostream>
#include <vector>
#include <chrono>
#include <string>
#include <functional>

// External Libraries
#include <H5Cpp.h>

// Local includes
#include "vec2.h"
#include "data_structures.h"
#include "helpers.h"
#include "integrators.h"
#include "seeds.h"

using namespace H5;

int main() {
    std::cout << "Starting line generation benchmark..." << std::endl;

    RawData data;
    readH5File("../data/isabel_2d.h5", &data);
    if (data.n_cols == 0 || data.n_rows == 0) {
        std::cerr << "Failed to load dataset or dataset is empty." << std::endl;
        return 1;
    }

    const int runs = 5;
    const float random_percentage = 0.01f;
    const size_t uniform_cell_distance = 10;
    const float even_line_distance = 2.0f;

    const float constant_step_size = 1.0f;
    const std::vector<size_t> max_steps_values = {100, 200, 400};

    const size_t constant_max_steps = 200;
    const std::vector<float> step_size_values = {0.25f, 0.5f, 1.0f};

    struct MethodProfile {
        std::string name;
        std::function<std::vector<Line>(float, size_t)> run;
    };

    std::vector<MethodProfile> methods = {
        {
            "random_rk4",
            [&](float step_size, size_t max_steps) {
                std::vector<Particle> seeds = getRandomSeed(&data, random_percentage);
                return rk4Integrator(&seeds, &data, step_size, max_steps);
            }
        },
        {
            "uniform_rk4",
            [&](float step_size, size_t max_steps) {
                std::vector<Particle> seeds = getUniformSeed(&data, uniform_cell_distance);
                return rk4Integrator(&seeds, &data, step_size, max_steps);
            }
        },
        {
            "even_seed",
            [&](float step_size, size_t max_steps) {
                return getEvenSeed(&data, even_line_distance, step_size, max_steps);
            }
        }
    };

    std::cout << "\nSweep 1: constant step_size=" << constant_step_size << ", vary max_steps" << std::endl;
    for (size_t max_steps : max_steps_values) {
        std::cout << "\nmax_steps=" << max_steps << std::endl;
        for (auto& method : methods) {
            long long sum_ms = 0;
            size_t last_line_count = 0;
            for (int run = 0; run < runs; ++run) {
                auto start = std::chrono::high_resolution_clock::now();
                std::vector<Line> lines = method.run(constant_step_size, max_steps);
                auto end = std::chrono::high_resolution_clock::now();
                sum_ms += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                last_line_count = lines.size();
            }
            std::cout << "  " << method.name << ": " << (sum_ms / runs) << " ms avg, " << last_line_count << " lines" << std::endl;
        }
    }

    std::cout << "\nSweep 2: constant max_steps=" << constant_max_steps << ", vary step_size" << std::endl;
    for (float step_size : step_size_values) {
        std::cout << "\nstep_size=" << step_size << std::endl;
        for (auto& method : methods) {
            long long sum_ms = 0;
            size_t last_line_count = 0;
            for (int run = 0; run < runs; ++run) {
                auto start = std::chrono::high_resolution_clock::now();
                std::vector<Line> lines = method.run(step_size, constant_max_steps);
                auto end = std::chrono::high_resolution_clock::now();
                sum_ms += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                last_line_count = lines.size();
            }
            std::cout << "  " << method.name << ": " << (sum_ms / runs) << " ms avg, " << last_line_count << " lines" << std::endl;
        }
    }

    return 0;
}
