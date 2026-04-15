#include "integrators.h"
#include "helpers.h"

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