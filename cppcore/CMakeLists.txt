cmake_minimum_required(VERSION 3.10)

# Set the project name
project(WannierHamiltonian VERSION 1.0 LANGUAGES CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Enable OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Set include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories("/home/jgarcia/projects/eigen-3.4.0/")

# Set the output directory for compiled objects
set(OUTPUT_DIR ${CMAKE_BINARY_DIR}/bin)
file(MAKE_DIRECTORY ${OUTPUT_DIR})

# Specify the source directory
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)

# Add the executable
add_executable(main ${SRC_DIR}/main.cpp)

# Link directories
link_directories(${CMAKE_SOURCE_DIR}/lib)

# Link libraries
target_link_libraries(main hamiltonian ${OpenMP_CXX_LIBRARIES})

# Specify the output directory for the executable
set_target_properties(main PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})

# Print a message after building the main program
add_custom_command(TARGET main POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E echo "Main program built successfully!")
