cmake_minimum_required (VERSION 2.6)
project (final C)

set(SOURCES twodheat.c)

include_directories(/usr/local/include)
add_executable(twodheat ${SOURCES})
target_link_libraries(twodheat umfpack amd cholmod colamd ccolamd camd metis
  suitesparseconfig cblas)