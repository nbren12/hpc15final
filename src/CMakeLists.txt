cmake_minimum_required (VERSION 2.6)
project (final)

add_library(laplace STATIC laplace.cpp)
target_link_libraries(laplace armadillo)


add_executable(test_heateq test_heateq.cpp)
target_link_libraries(test_heateq laplace)