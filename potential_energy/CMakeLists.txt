cmake_minimum_required(VERSION 3.15)
project(ar-h2-potential)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED  11)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra")

add_executable(
    main 
    main.cpp
)

find_package(GSL REQUIRED)
target_link_libraries(main GSL::gsl GSL::gslcblas)
