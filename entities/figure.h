#pragma once

#include <vector>
#include "triangle.h"
#include "vertex.h"
#include "material.h"

class Figure {
public:
    Figure() {
    }

    Figure(Material material, const std::vector<Vertex> &vertices, int index)
            : material_(material), index_(index) {
        int n = vertices.size();
        for (int i = 2; i < n; ++i) {
            triangles_.emplace_back(vertices[0], vertices[i - 1], vertices[i]);
        }
    }

    bool operator==(const Figure &rhs) const {
        return index_ == rhs.Index();
    }

    Material Mtl() const {
        return material_;
    }

    const std::vector<Triangle> &Triangles() const {
        return triangles_;
    }

    int Index() const {
        return index_;
    }

private:
    Material material_;
    std::vector<Triangle> triangles_;
    int index_;
};