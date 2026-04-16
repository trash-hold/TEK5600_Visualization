#include <iostream>
#include <vector>
#include <chrono>
#include <string>

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
    std::cout << "Starting integrator benchmarking..." << std::endl;

    RawData data;
    readH5File("../data/metsim1_2d.h5", &data);

    const std::vector<float> step_sizes = {0.1f, 1.0f};
    const float total_length = 20.0f;
    const int num_runs = 5;

    struct IntegratorInfo {
        const char* name;
        std::vector<Line> (*func)(const std::vector<Particle>*, const RawData*, const float&, const size_t&);
    } integrators[] = {
        {"euler", eulerIntegrator},
        {"rk4", rk4Integrator}
    };

    // Store total times for averaging
    std::vector<std::vector<long long>> total_times(step_sizes.size(), std::vector<long long>(2, 0)); // [step_size_index][integrator_index]

    for (int run = 0; run < num_runs; ++run) {
        std::cout << "Run " << (run + 1) << "/" << num_runs << std::endl;

        // Regenerate seeds for each run
        std::vector<Particle> seeds = getUniformSeed(&data, 5);

        for (size_t i = 0; i < step_sizes.size(); ++i) {
            float step_size = step_sizes[i];
            const float exact_steps = total_length / step_size;
            const size_t max_steps = static_cast<size_t>(std::round(exact_steps));

            if (std::fabs(exact_steps - static_cast<float>(max_steps)) > 1e-6f) {
                continue; // Skip invalid step sizes
            }

            for (size_t j = 0; j < 2; ++j) { // 0: euler, 1: rk4
                const auto& integrator = integrators[j];
                auto start = std::chrono::high_resolution_clock::now();
                std::vector<Line> field_lines = integrator.func(&seeds, &data, step_size, max_steps);
                auto end = std::chrono::high_resolution_clock::now();

                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
                total_times[i][j] += duration.count();
            }
        }
    }

    // Output averages
    std::cout << "\nAverage times over " << num_runs << " runs:" << std::endl;
    for (size_t i = 0; i < step_sizes.size(); ++i) {
        float step_size = step_sizes[i];
        const float exact_steps = total_length / step_size;
        const size_t max_steps = static_cast<size_t>(std::round(exact_steps));

        if (std::fabs(exact_steps - static_cast<float>(max_steps)) > 1e-6f) {
            continue;
        }

        std::cout << "step_size=" << step_size << " max_steps=" << max_steps << ":" << std::endl;
        for (size_t j = 0; j < 2; ++j) {
            const auto& integrator = integrators[j];
            long long avg_time = total_times[i][j] / num_runs;
            std::cout << "  " << integrator.name << ": " << avg_time << " ms" << std::endl;
        }
    }

    return 0;
}