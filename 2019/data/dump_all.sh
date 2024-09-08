#!/bin/bash

directory="./"
output_file="all_results.out"

if [ -f "$output_file" ]; then
  rm "$output_file"
fi

for file in "$directory"/*.out; do
  if [ -f "$file" ]; then
    # Append file name as a header
    echo "===== Dumping content from $file =====" >> "$output_file"
    
    cat "$file" >> "$output_file"
    
    echo -e "\n\n----- End of $file -----\n" >> "$output_file"
  fi
done

echo "All .out files have been concatenated into $output_file."
