#!/bin/bash

# Specify the folder path where your files are located
folder_path="./AAAI2020-benchmarks-cap1/"

# Convert all files to lowercase and apply replacements
# find "$folder_path" -type f -name '*.ltlf' -exec sh -c 'tr "[:upper:]" "[:lower:]" < "$1" > "$1.tmp" && mv "$1.tmp" "$1"' _ {} \;

# Replace 'g ' with 'G ' in all files in the folder
# find "$folder_path" -type f -name '*.ltlf' -exec sed -i '' 's/g /G /g' {} +

# Replace 'f ' with 'F ' in all files in the folder
# find "$folder_path" -type f -name '*.ltlf' -exec sed -i '' 's/f /F /g' {} +

# Replace 'x ' with 'X[!] ' in all files in the folder
# find "$folder_path" -type f -name '*.ltlf' -exec sed -i '' 's/x /X[!] /g' {} +

# find "$folder_path" -type f -name '*.ltlf' -exec sed -i '' 's/ 0/ temporaryvariable0/g' {} +

find "$folder_path" -type f -name '*.ltlf' -exec sed -i '' 's/ u / U /g' {} +

