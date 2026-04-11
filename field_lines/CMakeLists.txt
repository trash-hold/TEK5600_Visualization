cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 11)

project(ExerciseFieldLines)

find_package(HDF5 REQUIRED COMPONENTS CXX)

add_executable(ex1_test main.cpp)

target_include_directories(ex1_test PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(ex1_test PRIVATE ${HDF5_CXX_LIBRARIES})