#cmake -DCMAKE_BUILD_TYPE=Debug -S . -B build

cmake_minimum_required(VERSION 3.16)

project(particle CXX)
enable_testing()

set(CMAKE_CXX_STANDARD          17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)

add_executable(particle ParticleType.cpp ResonanceType.cpp Particle.cpp test2.cpp)
target_compile_options(particle PRIVATE -Wall -Wextra -g -fsanitize=address)
target_link_options(particle PRIVATE -fsanitize=address)