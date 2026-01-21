#!/bin/bash

# Start a timer for the entire script
echo "Starting cleanup process..."
time (

    # Remove the build directory
    echo "Removing 'build/' directory with verbose output..."
    time rm -rv build/
    echo "Removal finished."

)
echo "Total cleanup process finished."
