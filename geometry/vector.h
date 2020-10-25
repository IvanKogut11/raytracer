#pragma once

#include <cmath>

class Vector {
public:
    Vector() : x_(0.0), y_(0.0), z_(0.0) {
    }

    Vector(double x, double y, double z) : x_(x), y_(y), z_(z) {
    }

    bool operator==(const Vector& rhs) const {
        return x_ == rhs.x_ && y_ == rhs.y_ && z_ == rhs.z_;
    }

    double X() const {
        return x_;
    }

    double Y() const {
        return y_;
    }

    double Z() const {
        return z_;
    }

    double Len() const {
        return sqrt(x_ * x_ + y_ * y_ + z_ * z_);
    }

    Vector operator+(const Vector& rhs) const {
        return {x_ + rhs.X(), y_ + rhs.Y(), z_ + rhs.Z()};
    }

    Vector operator-(const Vector& rhs) const {
        return {x_ - rhs.X(), y_ - rhs.Y(), z_ - rhs.Z()};
    }

    Vector operator-() const {
        return {-x_, -y_, -z_};
    }

    Vector operator*(double c) const {
        return {c * x_, c * y_, c * z_};
    }

    Vector& operator+=(const Vector& rhs) {
        *this = *this + rhs;
        return *this;
    }

    Vector& operator-=(const Vector& rhs) {
        *this = *this - rhs;
        return *this;
    }

    Vector& operator*=(double c) {
        *this = *this * c;
        return *this;
    }

    Vector& Normalize() {
        double len = this->Len();
        x_ /= len;
        y_ /= len;
        z_ /= len;
        return *this;
    }

    Vector GetNormalized() const {
        auto copy = *this;
        return copy.Normalize();
    }

private:
    double x_;
    double y_;
    double z_;
};
