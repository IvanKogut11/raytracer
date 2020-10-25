cmake_minimum_required(VERSION 3.8)
project(raytracer)

set(CMAKE_CXX_STANDARD             17)
set(CMAKE_MODULE_PATH              "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_EXPORT_COMPILE_COMMANDS  ON)

find_package(PNG)
find_package(JPEG)

add_executable(raytracer main.cpp)

target_link_libraries(raytracer ${PNG_LIBRARY} ${JPEG_LIBRARIES})