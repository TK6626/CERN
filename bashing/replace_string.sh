#!/bin/bash

# Directory to search
DIR="/src/phi_phi_reconstruction"

# String to find and replace
SEARCH="/simple_analysis/"
REPLACE="/"

# Recursively replace in all files
find "$DIR" -type f -exec sed -i "s/${SEARCH}/${REPLACE}/g" {} +

echo "Replacement complete in all files under $DIR"
