#include "seeds.h"

std::vector<Particle> getRandomSeed(RawData* const data, float particle_percentage)
{
    const float eps = 1e-6;
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
        p.velocity = data->interpolate(Vec2(rand_col, rand_row));
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
        
        Vec2 velocity = data->interpolate(Vec2(rand_col, rand_row));

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
                        seed_queue.push_back(Particle{candidate, data->interpolate(candidate)});
                    }
                }
            }

            // Incrementing line ID for the next line
            lineID++;
        }
    }

    return field_lines;
}