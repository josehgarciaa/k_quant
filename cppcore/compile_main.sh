#!/bin/bash

# Set the compilation options (you can customize these as needed)
CXX="g++"                                # C++ compiler
CXX_FLAGS="-std=c++11 -Wall -fopenmp"    # Compiler flags with OpenMP support
INCLUDE_DIR="./include"                 # Path to the include directory
SRC_DIR="./src"                         # Path to the src directory
LIB_DIR="./lib"                         # Path to the lib directory
EIGEN_DIR="/usr/include/eigen3/"
OUTPUT_DIR="./bin"                      # Output directory for compiled objects

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Compile main.cpp
$CXX $CXX_FLAGS -I"$EIGEN_DIR" -I"$INCLUDE_DIR" -c main.cpp -o "$OUTPUT_DIR/main.o"

# Alternatively, if you used lib2.a, use this linking command instead:
$CXX -L"$LIB_DIR"  -Wl,-rpath,"$LIB_DIR" -fopenmp -o "$OUTPUT_DIR/main" "$OUTPUT_DIR/main.o" -lhamiltonian

echo "Main program built successfully! "
