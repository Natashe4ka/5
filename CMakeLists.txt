cmake_minimum_required(VERSION 2.8)

project(Ray_tracing)

set(Ray_tracing_src main.cpp vec3d.cpp create_image.cpp GetDataFromFile.cpp figure.cpp)
set(Ray_tracing_hdr head.hpp vec3d.hpp light.hpp figure.hpp figure_factory.hpp)

add_executable(${PROJECT_NAME} ${Ray_tracing_src} ${Ray_tracing_hdr})