#!/bin/bash

# Make a backup of the original source
cp src/reconstruct_all_data/rho_reconstruction/vary_cuts_glueball/dxy_glueball.cpp src/reconstruct_all_data/rho_reconstruction/vary_cuts_glueball/dxy_glueball.cpp.bak

for topology in "20" "40"
	do
	#loop through 20 and 40 topology 
	sed -i "s/TString topology = \".*\";/TString topology = \"$topology\";/"  src/reconstruct_all_data/rho_reconstruction/vary_cuts_glueball/dxy_glueball.cpp

	echo "Running with topology $topology" 
	for i in {0..4}
		do
	    # Replace the line with the new index
	    sed -i "s/Int_t index = .*/Int_t index = $i;/" src/reconstruct_all_data/rho_reconstruction/vary_cuts_glueball/dxy_glueball.cpp

	    # Compile
		g++ src/reconstruct_all_data/rho_reconstruction/vary_cuts_glueball/dxy_glueball.cpp lib/*.cpp -Llib -ldict `root-config --cflags --libs` -o bin/reconstruct_all_data/rho_reconstruction/vary_cuts_glueball/dxy_glueball 
	    # Run
	    echo "Running program with index=$i"
	    ./bin/reconstruct_all_data/rho_reconstruction/vary_cuts_glueball/dxy_glueball	
	done
done

# Restore original file
mv src/reconstruct_all_data/rho_reconstruction/vary_cuts_glueball/dxy_glueball.cpp.bak src/reconstruct_all_data/rho_reconstruction/vary_cuts_glueball/dxy_glueball.cpp
