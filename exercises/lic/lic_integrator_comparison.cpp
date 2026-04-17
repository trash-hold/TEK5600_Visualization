#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <filesystem>

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

using IntegratorFunc = std::vector<Line> (*)(const std::vector<Particle>*, const RawData*, const float&, const size_t&);

void runLICComparison(const RawData* data, const std::filesystem::path& output_dir, const std::string& dataset_name) {
    // Parameters for sweeps
    std::vector<float> step_sizes = {5.0f};
    size_t constant_kernel_length = 1;

    std::vector<size_t> kernel_lengths = {5, 10, 20};
    float constant_step_size = 5.0f;
    size_t num_runs = 5;

    // Integrators
    std::vector<std::pair<std::string, IntegratorFunc>> integrators = {
        {"euler", eulerIntegrator}
        //{"rk4", rk4Integrator}
    };

    // Sweep 1: Varying step_size, constant kernel_length
    std::cout << "Sweep 1: Varying step_size, kernel_length = " << constant_kernel_length << std::endl;
    for (float step_size : step_sizes) {
        for (auto& integrator_pair : integrators) {
            std::string integrator_name = integrator_pair.first;
            auto integrator = integrator_pair.second;

            // Measure time over num_runs
            double total_time = 0.0;
            for (size_t run = 0; run < num_runs; ++run) {
                auto start = std::chrono::high_resolution_clock::now();
                //PixelPlane result = simpleScaledLIC(data, 3, step_size, constant_kernel_length, integrator);
                PixelPlane result = simpleLIC(data, step_size, constant_kernel_length, integrator);
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                total_time += elapsed.count();
            }
            double avg_time = total_time / num_runs;
            std::cout << "  Step size: " << step_size << ", Integrator: " << integrator_name << ", Avg time: " << avg_time << " s" << std::endl;

            // Run once to save image
            PixelPlane result = simpleScaledLIC(data, 3, step_size, constant_kernel_length, integrator);
            sf::Image lic_image;
            createGSImage(&result, &lic_image);
            std::string filename = dataset_name + "_lic_step_" + std::to_string(step_size) + "_kl_" + std::to_string(constant_kernel_length) + "_" + integrator_name + ".png";
            if (!lic_image.saveToFile((output_dir / filename).string())) {
                std::cerr << "Failed to save image: " << filename << std::endl;
            }
        }
    }

    // Sweep 2: Varying kernel_length, constant step_size
    std::cout << "Sweep 2: Varying kernel_length, step_size = " << constant_step_size << std::endl;
    for (size_t kernel_length : kernel_lengths) {
        for (auto& integrator_pair : integrators) {
            std::string integrator_name = integrator_pair.first;
            auto integrator = integrator_pair.second;

            // Measure time over num_runs
            double total_time = 0.0;
            for (size_t run = 0; run < num_runs; ++run) {
                auto start = std::chrono::high_resolution_clock::now();
                //PixelPlane result = simpleScaledLIC(data, 3, constant_step_size, kernel_length, integrator);
                PixelPlane result = simpleLIC(data, constant_step_size, constant_kernel_length, integrator);
                
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                total_time += elapsed.count();
            }
            double avg_time = total_time / num_runs;
            std::cout << "  Kernel length: " << kernel_length << ", Integrator: " << integrator_name << ", Avg time: " << avg_time << " s" << std::endl;

            // Run once to save image
            PixelPlane result = simpleScaledLIC(data, 3, constant_step_size, kernel_length, integrator);
            sf::Image lic_image;
            createGSImage(&result, &lic_image);
            std::string filename = dataset_name + "_lic_step_" + std::to_string(constant_step_size) + "_kl_" + std::to_string(kernel_length) + "_" + integrator_name + ".png";
            if (!lic_image.saveToFile((output_dir / filename).string())) {
                std::cerr << "Failed to save image: " << filename << std::endl;
            }
        }
    }
}

int main() {
    std::cout << "Starting LIC Integrator Comparison..." << std::endl;

    RawData data;

    // Read file
    readH5File(isabel_file, &data);
    std::string dataset_name = "Kernel_Isabel";
    const std::filesystem::path output_dir = "../exercises/lic/img/kernel";
    std::filesystem::create_directories(output_dir);

    runLICComparison(&data, output_dir, dataset_name);

    std::cout << "Comparison complete. Images saved in " << output_dir.string() << std::endl;

    return 0;
}