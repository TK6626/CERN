#!/bin/bash

#use the base dir

base_dir=$1
file_name=$2
src_file="src/${base_dir}${file_name}.cpp" 

cp "$src_file" "${src_file}.bak"	
echo "running file $src_file"

for topology in "20" "40"
	do
	#loop through 20 and 40 topology 
	sed -i "s/TString topology = \".*\";/TString topology = \"$topology\";/" "$src_file"
	echo "Running Topology $topology"
	for i in {0..4}
		do
		sed -i "s/Int_t index = .*/Int_t index = $i;/" "$src_file"

	    # Compile and Run
	   	echo "Running program with index=$i"
		./bashing/run_file.sh "$base_dir$file_name"
	done
done
# Restore original file
mv "${src_file}.bak" "$src_file"


