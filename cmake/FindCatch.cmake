add_library(contrib_catch_main
  contrib/catch/catch_main.cpp ../geometry/vector.h ../entities/vertex.h ../parsing/obj_parser.h ../parsing/mtl_parser.h ../entities/material.h ../parsing/commons.h ../entities/sphere.h ../entities/light.h ../entities/triangle.h ../entities/figure.h ../geometry/geometry_calculations.h)

target_include_directories(contrib_catch_main
  PUBLIC contrib/catch)
