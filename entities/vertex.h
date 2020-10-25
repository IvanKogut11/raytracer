#pragma once

#include <cmath>
#include "../geometry/vector.h"

class Vertex {
public:
    bool is_normal_set;
    Vector vn;
    std::pair<double, double> vt;

    Vertex() {
    }

    Vertex(Vector point, bool is_normal_set = false, Vector vn = {},
           std::pair<double, double> vt = {})
            : is_normal_set(is_normal_set), vn(vn.GetNormalized()), vt(vt), point_(point) {
    }

    Vector Point() const {
        return point_;
    }

    bool IsNormalSet() {
        return is_normal_set;
    }

    Vector Normal() {
        return vn;
    }

    std::pair<double, double> VT() {
        return vt;
    }

private:
    Vector point_;
};