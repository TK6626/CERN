#!/bin/bash

# Directory to search
DIR="media/root_files/phi_phi_reconstruction/"

# Recursively remove all .bak files
find "$DIR" -type f -name "*.bak" -exec rm -f {} +

echo "All .bak files removed under $DIR"
