#pragma once

#include "../geometry/vector.h"

class Light {
public:
    Light(Vector point, Vector rgb) : point_(point), rgb_(rgb) {
    }

    Vector Point() const {
        return point_;
    }

    Vector RGB() const {
        return rgb_;
    }

private:
    Vector point_;
    Vector rgb_;
};