#!/bin/bash

#use base dir for cleaner writin
base_dir="reconstruct_all_data/rho_reconstruction/vary_cuts_glueball/" 



for file_name in "$@"
	do 
	src_file="src/${base_dir}${file_name}.cpp"
	
	# Make a backup of the original source
	cp "$src_file" "${src_file}.bak"
	
	echo "running file $src_file"

	for topology in "20" "40"
		do
		#loop through 20 and 40 topology 
		sed -i "s/TString topology = \".*\";/TString topology = \"$topology\";/" "$src_file"

		echo "Running with topology $topology" 
		for i in {0..4}
			do
	    	# Replace the line with the new index
	    	sed -i "s/Int_t index = .*/Int_t index = $i;/" "$src_file"

	    	# Compile and Run
	    	echo "Running program with index=$i"
			./bashing/run_file.sh "$base_dir$file_name"
		done
	done
done
# Restore original file
mv "${src_file}.bak" "$src_file"
