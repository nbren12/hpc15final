cmake_minimum_required (VERSION 2.6)
project (final)

include_directories(/usr/local/include)
add_library(heat SHARED laplace.cpp twodheat.cpp)
target_link_libraries(heat umfpack amd cholmod colamd ccolamd camd metis
  suitesparseconfig cblas)


add_executable(main main.cpp)
target_link_libraries(main heat)