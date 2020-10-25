#pragma once

#include "image.h"
#include "camera_options.h"
#include "render_options.h"
#include "parsing/obj_parser.h"
#include "entities/light.h"
#include "entities/sphere.h"
#include "entities/figure.h"
#include "geometry/geometry_calculations.h"

#include <string>
#include <vector>

std::vector<std::vector<Vector>> GetNotProcessedImage(const CameraOptions& camera_options,
                                                      const RenderOptions& render_options,
                                                      const std::vector<Figure>& figures,
                                                      const std::vector<Sphere>& spheres,
                                                      const std::vector<Light>& lights);

Image GetProcessedImage(const std::vector<std::vector<Vector>>& not_processed_image,
                        int image_width, int image_height);

Image Render(const std::string& filename, const CameraOptions& camera_options,
             const RenderOptions& render_options) {
    auto parser = ObjParser();
    auto figures = std::vector<Figure>();
    auto spheres = std::vector<Sphere>();
    auto lights = std::vector<Light>();
    parser.Parse(filename, &figures, &spheres, &lights);

    auto image = GetNotProcessedImage(camera_options, render_options, figures, spheres, lights);

    auto result_image =
        GetProcessedImage(image, camera_options.screen_width, camera_options.screen_height);

    return result_image;
}

Vector CalcDiffuse(const Vector& ray_to_light, const Vector& normal, const Vector& light_rgb) {
    return light_rgb *
           std::max(0., Geometry::Dot(normal.GetNormalized(), ray_to_light.GetNormalized()));
}

Vector CalcSpecular(const Vector& point, const Vector& ray_to_point, const Vector& ray_to_light,
                    const Vector& normal, const Vector& light_rgb, double ns_of_object) {
    Vector reflected_ray_to_light = Geometry::GetReflectedVector(-ray_to_light, normal);
    return light_rgb * pow(std::max(0., Geometry::Dot(-ray_to_point.GetNormalized(),
                                                      reflected_ray_to_light.GetNormalized())),
                           ns_of_object);
}

template <class Object>
Vector RecursiveRender(const Vector& point, Object object, const Material& mtl,
                       const Vector& ray_to_point, int current_depth,
                       const std::vector<Figure>& figures, const std::vector<Sphere>& spheres,
                       const std::vector<Light>& lights, int max_depth, bool is_in_object);

Vector GetReflectedIllumination(const Vector& point, const Vector& reflected_ray,
                                const Vector& normal_of_point, int current_depth,
                                const std::vector<Figure>& figures,
                                const std::vector<Sphere>& spheres,
                                const std::vector<Light>& lights, int max_depth, double eps) {
    Vector reflected_illumination;
    auto intersection_info = Geometry::TryGetRayAndClosestObjectIntersection(
        point + normal_of_point * eps, reflected_ray, figures, spheres);
    if (intersection_info.intersection_result != IntersectionResult::no_intersection) {
        if (intersection_info.intersection_result == IntersectionResult::with_sphere) {
            reflected_illumination = RecursiveRender(
                intersection_info.intersection_point, *intersection_info.intersected_sphere,
                intersection_info.intersected_sphere->Mtl(), reflected_ray, current_depth + 1,
                figures, spheres, lights, max_depth, false);
        } else {
            reflected_illumination = RecursiveRender(
                intersection_info.intersection_point, *intersection_info.intersected_triangle,
                intersection_info.intersected_figure->Mtl(), reflected_ray, current_depth + 1,
                figures, spheres, lights, max_depth, false);
        }
    }

    return reflected_illumination;
}

Vector GetRefractedIllumination(const Vector& point, const Vector& refracted_ray,
                                const Vector& normal_of_point, int current_depth,
                                const std::vector<Figure>& figures,
                                const std::vector<Sphere>& spheres,
                                const std::vector<Light>& lights, int max_depth, bool is_in_object,
                                double eps) {
    Vector refracted_illumination;
    auto intersection_info = Geometry::TryGetRayAndClosestObjectIntersection(
        point - normal_of_point * eps, refracted_ray, figures, spheres);
    if (intersection_info.intersection_result != IntersectionResult::no_intersection) {
        if (intersection_info.intersection_result == IntersectionResult::with_sphere) {
            refracted_illumination = RecursiveRender(
                intersection_info.intersection_point, *intersection_info.intersected_sphere,
                intersection_info.intersected_sphere->Mtl(), refracted_ray, current_depth + 1,
                figures, spheres, lights, max_depth, !is_in_object);
        } else {
            refracted_illumination = RecursiveRender(
                intersection_info.intersection_point, *intersection_info.intersected_triangle,
                intersection_info.intersected_figure->Mtl(), refracted_ray, current_depth + 1,
                figures, spheres, lights, max_depth, false);  // Maybe change
        }
    }
    return refracted_illumination;
}

template <class Object>
Vector RecursiveRender(const Vector& point, Object object, const Material& mtl,
                       const Vector& ray_to_point, int current_depth,
                       const std::vector<Figure>& figures, const std::vector<Sphere>& spheres,
                       const std::vector<Light>& lights, int max_depth, bool is_in_object) {
    const double eps = 1e-9;
    Vector illumination = mtl.Ka + mtl.Ke;
    auto normal_of_point = Geometry::FindNormalAtPoint(object, point).GetNormalized();
    if (!Geometry::AreNormalAndVectorToOneSide(normal_of_point, -ray_to_point.GetNormalized())) {
        normal_of_point = -normal_of_point;
    }
    Vector diffuse;
    Vector specular;
    for (const auto& light : lights) {
        auto ray_to_light = light.Point() - point;
        auto intersection_info = Geometry::TryGetRayAndClosestObjectIntersection(
            point + normal_of_point * eps, ray_to_light.GetNormalized(), figures, spheres);
        if (intersection_info.intersection_result != IntersectionResult::no_intersection &&
            (intersection_info.intersection_point - point).Len() < ray_to_light.Len()) {
            continue;
        }
        diffuse += CalcDiffuse(ray_to_light, normal_of_point, light.RGB());
        specular +=
            CalcSpecular(point, ray_to_point, ray_to_light, normal_of_point, light.RGB(), mtl.Ns);
    }

    if (current_depth != max_depth && mtl.illum > 2) {
        if (!is_in_object) {
            Vector reflected_ray = Geometry::GetReflectedVector(ray_to_point, normal_of_point);
            auto reflected_illumination =
                GetReflectedIllumination(point, reflected_ray, normal_of_point, current_depth,
                                         figures, spheres, lights, max_depth, eps);
            diffuse += CalcDiffuse(reflected_ray, normal_of_point, reflected_illumination);
            specular += CalcSpecular(point, ray_to_point, reflected_ray, normal_of_point,
                                     reflected_illumination, mtl.Ns);
        }

        double n_out = 1, n_in = mtl.Ni;
        if (is_in_object) {
            std::swap(n_out, n_in);
        }
        Vector refracted_ray =
            Geometry::GetRefractedVector(ray_to_point, n_out, n_in, normal_of_point);
        auto refracted_illumination =
            GetRefractedIllumination(point, refracted_ray, normal_of_point, current_depth, figures,
                                     spheres, lights, max_depth, is_in_object, eps);
        if (is_in_object) {
            illumination += refracted_illumination;
        } else {
            illumination += refracted_illumination * (1.0 - mtl.d);
        }
    }

    illumination += {mtl.Kd.X() * diffuse.X() + mtl.Ks.X() * specular.X(),
                     mtl.Kd.Y() * diffuse.Y() + mtl.Ks.Y() * specular.Y(),
                     mtl.Kd.Z() * diffuse.Z() + mtl.Ks.Z() * specular.Z()};
    return illumination;
}

std::vector<std::vector<Vector>> GetNotProcessedImage(const CameraOptions& camera_options,
                                                      const RenderOptions& render_options,
                                                      const std::vector<Figure>& figures,
                                                      const std::vector<Sphere>& spheres,
                                                      const std::vector<Light>& lights) {
    double scale = tan(camera_options.fov * 0.5);
    double image_aspect_ratio =
        camera_options.screen_width / static_cast<double>(camera_options.screen_height);
    auto result_image = std::vector<std::vector<Vector>>(
        camera_options.screen_height, std::vector<Vector>(camera_options.screen_width, {0, 0, 0}));

    auto look_from = Vector(camera_options.look_from[0], camera_options.look_from[1],
                            camera_options.look_from[2]);
    auto look_to =
        Vector(camera_options.look_to[0], camera_options.look_to[1], camera_options.look_to[2]);
    for (int i = 0; i < camera_options.screen_height; ++i) {
        for (int j = 0; j < camera_options.screen_width; ++j) {
            double x = (2 * (j + 0.5) / static_cast<double>(camera_options.screen_width) - 1) *
                       image_aspect_ratio * scale;
            double y =
                (1 - 2 * (i + 0.5) / static_cast<double>(camera_options.screen_height)) * scale;
            auto point = Vector(x, y, -1);
            auto ray_in_space =
                Geometry::GetVectorFromCameraToSpace(look_from, look_to, point.GetNormalized())
                    .GetNormalized();

            IntersectionInfo intersection_info = Geometry::TryGetRayAndClosestObjectIntersection(
                look_from, ray_in_space, figures, spheres);
            if (intersection_info.intersection_result == IntersectionResult::with_figure) {
                result_image[i][j] = RecursiveRender(
                    intersection_info.intersection_point, *intersection_info.intersected_triangle,
                    intersection_info.intersected_figure->Mtl(), ray_in_space, 1, figures, spheres,
                    lights, render_options.depth, false);
            } else if (intersection_info.intersection_result == IntersectionResult::with_sphere) {
                result_image[i][j] = RecursiveRender(
                    intersection_info.intersection_point, *intersection_info.intersected_sphere,
                    intersection_info.intersected_sphere->Mtl(), ray_in_space, 1, figures, spheres,
                    lights, render_options.depth, false);
            }
        }
    }

    return result_image;
}

void DoToneMapping(std::vector<std::vector<Vector>>* not_processed_image, int image_width,
                   int image_height) {
    double max_c = 0;
    for (int i = 0; i < image_height; ++i) {
        for (int j = 0; j < image_width; ++j) {
            max_c = std::max(max_c, (*not_processed_image)[i][j].X());
            max_c = std::max(max_c, (*not_processed_image)[i][j].Y());
            max_c = std::max(max_c, (*not_processed_image)[i][j].Z());
        }
    }

    for (int i = 0; i < image_height; ++i) {
        for (int j = 0; j < image_width; ++j) {
            double r = (*not_processed_image)[i][j].X();
            double g = (*not_processed_image)[i][j].Y();
            double b = (*not_processed_image)[i][j].Z();
            r = (r * (1 + r / (max_c * max_c))) / (1 + r);
            g = (g * (1 + g / (max_c * max_c))) / (1 + g);
            b = (b * (1 + b / (max_c * max_c))) / (1 + b);
            (*not_processed_image)[i][j] = Vector(r, g, b);
        }
    }
}

void DoGammaCorrection(std::vector<std::vector<Vector>>* not_processed_image, int image_width,
                       int image_height) {
    const double degree = 1. / 2.2;
    for (int i = 0; i < image_height; ++i) {
        for (int j = 0; j < image_width; ++j) {
            double r = pow((*not_processed_image)[i][j].X(), degree);
            double g = pow((*not_processed_image)[i][j].Y(), degree);
            double b = pow((*not_processed_image)[i][j].Z(), degree);
            (*not_processed_image)[i][j] = Vector(r, g, b);
        }
    }
}

Image GetProcessedImage(const std::vector<std::vector<Vector>>& not_processed_image,
                        int image_width, int image_height) {
    auto copy_not_processed_image = not_processed_image;
    auto image = Image(image_width, image_height);
    DoToneMapping(&copy_not_processed_image, image_width, image_height);
    DoGammaCorrection(&copy_not_processed_image, image_width, image_height);
    for (int i = 0; i < image_height; ++i) {
        for (int j = 0; j < image_width; ++j) {
            RGB rgb;
            rgb.r = copy_not_processed_image[i][j].X() * 255;
            rgb.g = copy_not_processed_image[i][j].Y() * 255;
            rgb.b = copy_not_processed_image[i][j].Z() * 255;
            image.SetPixel(rgb, i, j);
        }
    }
    return image;
}