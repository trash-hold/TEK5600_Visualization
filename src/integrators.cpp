#include "integrators.h"
#include "helpers.h"

#define EPS_SQ 1e-8f
#define EPS_POS 1e-4f

Particle euler_step(const Vec2 &pos, const RawData *data, const float &step_size)
{
    Vec2 velocity = data->interpolate(pos).normalized() * step_size;
    Vec2 new_pos = pos + velocity;
    return Particle{new_pos, data->interpolate(new_pos)};
}

std::vector<Line> eulerIntegrator(const std::vector<Particle> *seeds, const RawData *data, const float &step_size, const size_t &max_steps)
{
    std::vector<Line> field_lines;
    field_lines.reserve(seeds->size()); // We compute field-lines for all generated seeds

    size_t n_rows = data->n_rows;
    size_t n_cols = data->n_cols;

    for (const auto &seed : *seeds)
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

        for (int i = 0; i < max_steps && (f_inbounds || b_inbounds); i++)
        {

            if (f_inbounds)
            { // To avoid checking overhead, the f_inbounds is updated here
                if (data->isValid(static_cast<size_t>(f_current_pos.y), static_cast<size_t>(f_current_pos.x)) == true)
                {
                    Particle forward = euler_step(f_current_pos, data, step_size);
                    f_current_pos = forward.position;
                    // Add the new point to the line
                    deq_line.push_back(forward);

                    // Checking if the field isn't to weak
                    if (forward.velocity.x * forward.velocity.x + forward.velocity.y * forward.velocity.y < EPS_SQ)
                        f_inbounds = false;
                }
                else
                    f_inbounds = false;
            }

            if (b_inbounds)
            {
                if (data->isValid(static_cast<size_t>(b_current_pos.y), static_cast<size_t>(b_current_pos.x)))
                {
                    Particle backward = euler_step(b_current_pos, data, -step_size);
                    b_current_pos = backward.position;
                    // Add the new point to the line
                    deq_line.push_front(backward);

                    // Checking if the field isn't to weak
                    if (backward.velocity.x * backward.velocity.x + backward.velocity.y * backward.velocity.y < EPS_SQ)
                        b_inbounds = false;
                }
                else
                    b_inbounds = false;
            }
        }

        // Deque to vector
        std::vector<Particle> particles(deq_line.begin(), deq_line.end());
        field_lines.push_back(Line{std::move(particles)});
    }
    return field_lines;
}

Particle rk4_step(const Vec2 &pos, const RawData *data, const float &step_size)
{
    Vec2 k1 = data->interpolate(pos);
    // Premature check if the field isn't too weak
    if(k1.x*k1.x + k1.y*k1.y < EPS_SQ)
        return Particle{pos, k1};

    k1 = k1.normalized() * step_size;
    Vec2 k2 = data->interpolate(pos + k1 * 0.5f).normalized() * step_size;
    Vec2 k3 = data->interpolate(pos + k2 * 0.5f).normalized() * step_size;
    Vec2 k4 = data->interpolate(pos + k3).normalized() * step_size;

    Vec2 new_pos = pos + (k1 + k2 * 2.0f + k3 * 2.0f + k4) * (1.0f / 6.0f);
    return Particle{new_pos, data->interpolate(new_pos)};
}

std::vector<Line> rk4Integrator(const std::vector<Particle> *seeds, const RawData *data, const float &step_size, const size_t &max_steps)
{
    std::vector<Line> field_lines;
    field_lines.reserve(seeds->size()); // We compute field-lines for all generated seeds

    size_t n_rows = data->n_rows;
    size_t n_cols = data->n_cols;

    for (const auto &seed : *seeds)
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

        for (int i = 0; i < max_steps && (f_inbounds || b_inbounds); i++)
        {
            if (f_inbounds)
            { // To avoid checking overhead, the f_inbounds is updated here
                if (data->isValid(static_cast<size_t>(f_current_pos.y), static_cast<size_t>(f_current_pos.x)) == true)
                {
                    Particle forward = rk4_step(f_current_pos, data, step_size);
                    f_current_pos = forward.position;

                    // Check for weak field
                    if (forward.velocity.x * forward.velocity.x + forward.velocity.y*forward.velocity.y < EPS_SQ)
                        f_inbounds = false;

                    else
                        deq_line.push_back(forward);

                }
                else
                    f_inbounds = false;
            }

            if (b_inbounds)
            {
                if (data->isValid(static_cast<size_t>(b_current_pos.y), static_cast<size_t>(b_current_pos.x)))
                {
                    Particle backward = rk4_step(b_current_pos, data, -step_size);
                    b_current_pos = backward.position;

                    // Check for weak field
                    if (backward.velocity.x * backward.velocity.x + backward.velocity.y * backward.velocity.y < EPS_SQ)
                        b_inbounds = false;
                    else
                        deq_line.push_front(backward); // Add the new point to the line
                }
                else
                    b_inbounds = false;
            }
        }

        // Deque to vector
        std::vector<Particle> particles(deq_line.begin(), deq_line.end());
        field_lines.push_back(Line{std::move(particles)});
    }
    return field_lines;
}