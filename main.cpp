#include <iostream>
#include "raytracer.h"

int main(int args_cnt, char** args) {
    auto obj_filename = std::string(args[1]);

    auto camera_width = std::atoi(args[2]);
    auto camera_height = std::atoi(args[3]);
    CameraOptions camera_opts(camera_width, camera_height);
    double from_x = std::atof(args[4]), from_y = std::atof(args[5]), from_z = std::atof(args[6]);
    double to_x = std::atof(args[7]), to_y = std::atof(args[8]), to_z = std::atof(args[9]);
    camera_opts.look_from = std::array<double, 3>{from_x, from_y, from_z};
    camera_opts.look_to = std::array<double, 3>{to_x, to_y, to_z};
    RenderOptions render_opts{std::atoi(args[10])};

    auto rendered_filename = std::string(args[11]);
    auto rendered = Render(obj_filename, camera_opts, render_opts);
    rendered.Write(rendered_filename);
    return 0;
}
