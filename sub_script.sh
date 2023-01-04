#!/bin/bash

# Initialize and Load Modules
# no modules to load when using MATLAB

echo "My task ID: " $LLSUB_RANK
echo "Number of Tasks: " $LLSUB_SIZE

matlab -batch "MyTaskID=$LLSUB_RANK;NumberOfTasks=$LLSUB_SIZE;lin_global_search_supercloud"