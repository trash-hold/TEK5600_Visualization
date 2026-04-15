cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 17)

project(ExerciseFieldLines)

find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(SFML 3 REQUIRED COMPONENTS Graphics)

# Get all subdirectories in exercises/
file(GLOB EXERCISE_DIRS LIST_DIRECTORIES true "exercises/*")

# For each exercise directory, create an executable
foreach(EXERCISE_DIR ${EXERCISE_DIRS})
    if(IS_DIRECTORY ${EXERCISE_DIR})
        get_filename_component(EXERCISE_NAME ${EXERCISE_DIR} NAME)
        
        # Find the main.cpp file in the exercise directory
        file(GLOB MAIN_FILE "${EXERCISE_DIR}/main.cpp")
        
        if(MAIN_FILE)
            add_executable(${EXERCISE_NAME}_exercise 
                ${MAIN_FILE}
                src/helpers.cpp
                src/integrators.cpp
                src/seeds.cpp
                src/lic.cpp
                src/visualization.cpp
            )
            
            target_include_directories(${EXERCISE_NAME}_exercise PRIVATE ${HDF5_INCLUDE_DIRS} inc/)
            target_link_libraries(${EXERCISE_NAME}_exercise PRIVATE ${HDF5_CXX_LIBRARIES} SFML::Graphics)
        endif()
    endif()
endforeach()

# Original combined exercise (optional)
add_executable(ex1_test 
    main.cpp
    src/helpers.cpp
    src/integrators.cpp
    src/seeds.cpp
    src/lic.cpp
    src/visualization.cpp
)

target_include_directories(ex1_test PRIVATE ${HDF5_INCLUDE_DIRS} inc/)
target_link_libraries(ex1_test PRIVATE ${HDF5_CXX_LIBRARIES} SFML::Graphics)