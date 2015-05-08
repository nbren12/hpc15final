cmake_minimum_required (VERSION 2.6)
project (final)

add_executable(main twodheat.cpp laplace.cpp)
target_link_libraries(main armadillo)