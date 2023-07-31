#!/bin/bash

# Set the compilation options (you can customize these as needed)
CXX="g++"                          # C++ compiler
CXX_FLAGS="-std=c++11 -Wall -O3 -fopenmp"       # Compiler flags
INCLUDE_DIR="./include"           # Path to the include directory
SRC_DIR="./src"                   # Path to the src directory
OUTPUT_DIR="./bin"                # Output directory for compiled objects
LIB_DIR="./lib"                   # Output directory for compiled libraries
rm $LIB_DIR/*
# Create the output directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LIB_DIR"

$CXX $CXX_FLAGS -I"$INCLUDE_DIR" -c "$SRC_DIR/sparse_matrices.cpp" -o "$OUTPUT_DIR/libsparse_matrices.o"
$CXX $CXX_FLAGS -I"$INCLUDE_DIR" -c "$SRC_DIR/wannier90_utils.cpp" -o "$OUTPUT_DIR/libwannier90_utils.o"
$CXX $CXX_FLAGS -I"$INCLUDE_DIR" -c "$SRC_DIR/hamiltonian.cpp" -o "$OUTPUT_DIR/libhamiltonian.o"

# Create a static library:
ar rcs "$LIB_DIR/libhamiltonian.a" "$OUTPUT_DIR/libhamiltonian.o" "$OUTPUT_DIR/libsparse_matrices.o" "$OUTPUT_DIR/libwannier90_utils.o"

# Clean up intermediate object files
rm -rf "$OUTPUT_DIR"

echo "Libraries compiled successfully!"
