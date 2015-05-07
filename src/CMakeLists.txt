cmake_minimum_required (VERSION 2.6)
project (final)

set(SOURCES laplace.cpp)

include_directories(/usr/local/include)
add_library(twodheat SHARED ${SOURCES})
target_link_libraries(twodheat umfpack amd cholmod colamd ccolamd camd metis
  suitesparseconfig cblas)


add_executable(main main.cpp)
target_link_libraries(main twodheat)