#pragma once

#include <cmath>
#include <vector>
#include "../entities/triangle.h"
#include "../entities/sphere.h"
#include "vector.h"

enum class IntersectionResult {
    no_intersection, with_sphere, with_figure
};

struct IntersectionInfo {
    IntersectionResult intersection_result = IntersectionResult::no_intersection;
    Vector intersection_point;
    const Sphere *intersected_sphere{};
    const Figure *intersected_figure{};
    const Triangle *intersected_triangle{};

    IntersectionInfo() = default;
};

class Geometry {
public:
    static double Dot(Vector v1, Vector v2) {
        return v1.X() * v2.X() + v1.Y() * v2.Y() + v1.Z() * v2.Z();
    }

    static Vector Cross(Vector v1, Vector v2) {
        return {v1.Y() * v2.Z() - v1.Z() * v2.Y(), -(v1.X() * v2.Z() - v1.Z() * v2.X()),
                v1.X() * v2.Y() - v1.Y() * v2.X()};
    }

    static double Triple(Vector v1, Vector v2, Vector v3) {
        return Dot(v1, Cross(v2, v3));
    }

    static bool IsPointInTriangle(const Triangle &tr, const Vector &pt) {
        auto n =
                Cross(tr[1].Point() - tr[0].Point(), tr[2].Point() - tr[0].Point()).GetNormalized();
        auto dot1 = Dot(n, Cross(tr[1].Point() - tr[0].Point(), pt - tr[0].Point()));
        auto dot2 = Dot(n, Cross(tr[2].Point() - tr[1].Point(), pt - tr[1].Point()));
        auto dot3 = Dot(n, Cross(tr[0].Point() - tr[2].Point(), pt - tr[2].Point()));
        return DoubleLessEqual(dot1, 0) && DoubleLessEqual(dot2, 0) && DoubleLessEqual(dot3, 0) ||
               DoubleLessEqual(0, dot1) && DoubleLessEqual(0, dot2) && DoubleLessEqual(0, dot3);
    }

    static bool TryFindObjectAndRayIntersection(const Triangle &tr, const Vector &st,
                                                const Vector &dir, Vector *inter) {
        auto n = Cross(tr[1].Point() - tr[0].Point(), tr[2].Point() - tr[0].Point());
        if (!AreNormalAndVectorToOneSide(n, -dir)) {
            n = -n;
        }
        auto dot1 = Dot(dir, n);
        if (DoubleEqual(dot1, 0)) {
            return false;
        }
        auto dot2 = Dot(st - tr[0].Point(), n);
        auto t = -dot2 / dot1;
        auto p_in_plane = st + dir * t;
        if (DoubleGreat(t, 0) && IsPointInTriangle(tr, p_in_plane)) {
            (*inter) = p_in_plane;
            return true;
        }
        return false;
    }

    static bool TryFindObjectAndRayIntersection(const Sphere &sphere, const Vector &st,
                                                const Vector &dir, Vector *inter) {
        std::pair<double, double> roots;
        auto a = Dot(dir, dir);
        auto b = 2 * Dot(dir, st - sphere.Center());
        auto c =
                Dot(st - sphere.Center(), st - sphere.Center()) - sphere.Radius() * sphere.Radius();
        if (!TrySolveQuadraticEquation(a, b, c, &roots)) {
            return false;
        }
        if (DoubleGreat(roots.first, 0)) {
            (*inter) = st + dir * roots.first;
            return true;
        }
        if (DoubleGreat(roots.second, 0)) {
            (*inter) = st + dir * roots.second;
            return true;
        }
        return false;
    }

    static double FindTriangleSquare(const Triangle &tr) {
        return 0.5 * Cross(tr[1].Point() - tr[0].Point(), tr[2].Point() - tr[0].Point()).Len();
    }

    static Vector FindNormalAtPoint(const Sphere &sphere, const Vector &point) {
        return point - sphere.Center();
    }

    static Vector FindNormalAtPoint(const Triangle &tr, const Vector &point) {
        if (!tr[0].is_normal_set || !tr[1].is_normal_set || !tr[2].is_normal_set) {
            return Cross(tr[1].Point() - tr[0].Point(), tr[2].Point() - tr[0].Point());
        }
        auto tr_square = FindTriangleSquare(tr);
        auto square_0_1 = FindTriangleSquare({tr[0].Point(), tr[1].Point(), point});
        auto square_0_2 = FindTriangleSquare({tr[0].Point(), tr[2].Point(), point});
        auto square_1_2 = FindTriangleSquare({tr[1].Point(), tr[2].Point(), point});
        return tr[0].Normal() * (square_1_2 / tr_square) +
               tr[1].Normal() * (square_0_2 / tr_square) +
               tr[2].Normal() * (square_0_1 / tr_square);
    }

    // From Link
    static Vector GetVectorFromCameraToSpace(const Vector &from, const Vector &to,
                                             const Vector &v) {
        auto matrix = GetLookAtMatrix(from, to);
        auto x = v.X() * matrix[0][0] + v.Y() * matrix[1][0] + v.Z() * matrix[2][0];
        auto y = v.X() * matrix[0][1] + v.Y() * matrix[1][1] + v.Z() * matrix[2][1];
        auto z = v.X() * matrix[0][2] + v.Y() * matrix[1][2] + v.Z() * matrix[2][2];

        return Vector(x, y, z);
    }

    static bool AreNormalAndVectorToOneSide(const Vector &normal, const Vector &vector) {
        return DoubleLessEqual(0, Dot(normal, vector));
    }

    static IntersectionInfo TryGetRayAndClosestObjectIntersection(
            const Vector &point, const Vector &ray, const std::vector<Figure> &figures,
            const std::vector<Sphere> &spheres) {
        double min_dist;
        IntersectionInfo intersection_info;
        intersection_info.intersection_result = IntersectionResult::no_intersection;

        for (auto &figure : figures) {
            for (auto &tr : figure.Triangles()) {
                Vector intersection_point;
                if (TryFindObjectAndRayIntersection(tr, point, ray.GetNormalized(),
                                                    &intersection_point)) {
                    if (intersection_info.intersection_result ==
                        IntersectionResult::no_intersection ||
                        DoubleGreat(min_dist, (intersection_point - point).Len())) {
                        min_dist = (intersection_point - point).Len();
                        intersection_info.intersection_point = intersection_point;
                        intersection_info.intersection_result = IntersectionResult::with_figure;
                        intersection_info.intersected_figure = &figure;
                        intersection_info.intersected_triangle = &tr;
                    }
                }
            }
        }

        for (auto &sphere : spheres) {
            Vector intersection_point;
            if (TryFindObjectAndRayIntersection(sphere, point, ray.GetNormalized(),
                                                &intersection_point)) {
                if (intersection_info.intersection_result == IntersectionResult::no_intersection ||
                    DoubleGreat(min_dist, (intersection_point - point).Len())) {
                    min_dist = (intersection_point - point).Len();
                    intersection_info.intersection_point = intersection_point;
                    intersection_info.intersection_result = IntersectionResult::with_sphere;
                    intersection_info.intersected_sphere = &sphere;
                }
            }
        }

        return intersection_info;
    }

    static Vector GetReflectedVector(const Vector &vector, const Vector &normal) {
        auto n_normalized = normal.GetNormalized();
        return n_normalized * (-Dot(vector, n_normalized) * 2) + vector;
    }

    static Vector GetRefractedVector(const Vector &vector, double n_out, double n_in,
                                     const Vector &normal) {
        auto v_ray = vector.GetNormalized() * n_out;
        auto n_normalized = normal.GetNormalized();
        auto coef = sqrt((n_in * n_in - n_out * n_out) /
                         (Dot(v_ray, n_normalized) * Dot(v_ray, n_normalized)) +
                         1) -
                    1;
        return v_ray + n_normalized * Dot(v_ray, n_normalized) * coef;
    }

private:
    static double DoubleEqual(double a, double b) {
        const double eps = 1e-9;
        return fabs(a - b) < eps;
    }

    static double DoubleLessEqual(double a, double b) {
        return a < b || DoubleEqual(a, b);
    }

    static double DoubleGreat(double a, double b) {
        return a > b && !DoubleEqual(a, b);
    }

    static bool TrySolveQuadraticEquation(double a, double b, double c,
                                          std::pair<double, double> *roots) {
        auto discriminant = b * b - 4 * a * c;
        if (DoubleGreat(0, discriminant)) {
            return false;
        }
        auto root1 = (-b + sqrt(discriminant)) / (2 * a);
        auto root2 = (-b - sqrt(discriminant)) / (2 * a);
        if (DoubleGreat(root2, root1)) {
            (*roots) = {root1, root2};
        } else {
            (*roots) = {root2, root1};
        }

        return true;
    }

    // From Link
    static std::vector<std::vector<double>> GetLookAtMatrix(const Vector &from, const Vector &to,
                                                            const Vector &tmp = Vector(0, 1, 0)) {
        Vector forward = (from - to).GetNormalized();
        Vector right = Cross(tmp.GetNormalized(), forward);
        Vector up = Cross(forward, right);
        if (forward == Vector(0, 1, 0)) {
            right = Vector(1, 0, 0);
            up = Vector(0, 0, -1);
        } else if (forward == Vector(0, -1, 0)) {
            right = Vector(-1, 0, 0);
            up = Vector(0, 0, 1);
        }
        forward.Normalize();
        right.Normalize();
        up.Normalize();
        auto matrix = std::vector<std::vector<double>>(4, std::vector<double>(4));
        matrix[0][0] = right.X();
        matrix[0][1] = right.Y();
        matrix[0][2] = right.Z();
        matrix[1][0] = up.X();
        matrix[1][1] = up.Y();
        matrix[1][2] = up.Z();
        matrix[2][0] = forward.X();
        matrix[2][1] = forward.Y();
        matrix[2][2] = forward.Z();

        matrix[3][0] = from.X();
        matrix[3][1] = from.Y();
        matrix[3][2] = from.Z();

        return matrix;
    }
};