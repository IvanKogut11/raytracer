#pragma once

#include "../geometry/vector.h"
#include "material.h"
#include "figure.h"

class Sphere {
public:
    Sphere() {
    }

    Sphere(Material material, Vector center, double radius, int index)
            : material_(material), center_(center), radius_(radius), index_(index) {
    }

    bool operator==(const Sphere &rhs) const {
        return index_ == rhs.Index();
    }

    bool operator==(const Figure &rhs) const {
        return false;
    }

    bool operator==(const Triangle &rhs) const {
        return false;
    }

    Material Mtl() const {
        return material_;
    }

    Vector Center() const {
        return center_;
    }

    double Radius() const {
        return radius_;
    }

    int Index() const {
        return index_;
    }

private:
    Material material_;
    Vector center_;
    double radius_;
    int index_;
};