#pragma once

#include "vertex.h"
#include <exception>

class Triangle {
public:
    Triangle() {
    }

    Triangle(Vertex p1, Vertex p2, Vertex p3) : p1_(p1), p2_(p2), p3_(p3) {
    }

    Vertex operator[](int index) const {
        if (index == 0) {
            return p1_;
        } else if (index == 1) {
            return p2_;
        } else if (index == 2) {
            return p3_;
        }
        throw std::exception();
    }

private:
    Vertex p1_;
    Vertex p2_;
    Vertex p3_;
};