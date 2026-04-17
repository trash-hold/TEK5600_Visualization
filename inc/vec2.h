#pragma once

#include <cmath>

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