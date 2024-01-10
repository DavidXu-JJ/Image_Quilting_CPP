cmake_minimum_required(VERSION 3.20)
project(image_quilting VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(image_quilting main.cpp)

target_include_directories(image_quilting PRIVATE third_party third_party/eigen)