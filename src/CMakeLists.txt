cmake_minimum_required (VERSION 2.6)
project (final C)

set(SOURCES main.c)

include_directories(/usr/local/include)
add_executable(main ${SOURCES})