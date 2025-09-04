#!/bin/bash

# take the argument to compile and run the file 
src_file="$1"
src_path="src/${src_file}.cpp"
bin_path="bin/${src_file}"

# compile file with all the helper files in the library
g++ "$src_path" lib/*.cpp -Llib -ldict `root-config --cflags --libs` -lRooFit -lRooFitCore -o "$bin_path"

#run the binary
"$bin_path"
