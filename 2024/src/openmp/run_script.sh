#!/bin/bash

# Set the output file
OUTPUT_FILE="./logs/results_$(date +%Y%m%d%H%M%S).log"
MESSAGE=$1

# Run the program multiple times, appending results to the file
echo "---*-----*------*------*--------*-----*--------" >> $OUTPUT_FILE
echo $MESSAGE >> $OUTPUT_FILE
for i in {1..5}; do
    echo "Run #$i" >> $OUTPUT_FILE   # Add a message indicating the run number
    ../../bin/openmp/main_openmp >> $OUTPUT_FILE   # Run the program and append the output
    echo "Completed run #$i" >> $OUTPUT_FILE  # Optional: Add a completion message
    echo "---------------------------------" >> $OUTPUT_FILE  # Separate the runs with a line
done

echo "All runs completed. Results are stored in $OUTPUT_FILE."
