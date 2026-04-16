cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 17)

project(ExerciseFieldLines)

find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(SFML 3 REQUIRED COMPONENTS Graphics)

# Get all subdirectories in exercises/
file(GLOB EXERCISE_DIRS LIST_DIRECTORIES true "exercises/*")

# For each exercise directory, create an executable for every .cpp file in that directory
foreach(EXERCISE_DIR ${EXERCISE_DIRS})
    if(IS_DIRECTORY ${EXERCISE_DIR})
        get_filename_component(EXERCISE_NAME ${EXERCISE_DIR} NAME)

        file(GLOB CPP_FILES "${EXERCISE_DIR}/*.cpp")
        foreach(CPP_FILE ${CPP_FILES})
            get_filename_component(CPP_NAME ${CPP_FILE} NAME_WE)
            set(EXE_NAME "${EXERCISE_NAME}_${CPP_NAME}")

            add_executable(${EXE_NAME}
                ${CPP_FILE}
                src/helpers.cpp
                src/integrators.cpp
                src/seeds.cpp
                src/lic.cpp
                src/visualization.cpp
            )

            target_include_directories(${EXE_NAME} PRIVATE ${HDF5_INCLUDE_DIRS} inc/)
            target_link_libraries(${EXE_NAME} PRIVATE ${HDF5_CXX_LIBRARIES} SFML::Graphics)
        endforeach()
    endif()
endforeach()