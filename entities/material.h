#pragma once

#include "../geometry/vector.h"

struct Material {
    double illum = 0;
    Vector Ka = {0., 0., 0.}, Kd = {0., 0., 0.}, Ks = {0., 0., 0.}, Ke = {0., 0., 0.};
    double Ns = 0, Ni = 0, d, Tr;

    bool IsNotTransparent() {
        return d == 1 || Tr == 0;
    }
};