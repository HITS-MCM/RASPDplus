#!/usr/bin/env bash

# Fill both empty variables with the correct values for RASPD+ to work
# Alternatively you can place the two variable exports in your ~/.bashrc or ~/.bashrc_profile
export raspd_root='' # Root of the cloned RASPD+ git repository
export conda_root='' # Path to your conda installation eg, /home/my_user/miniconda3
export TRAPP=''  # Path of the TRAPP installation (download TRAPP from https://www.h-its.org/downloads/trapp/)
# for verbosity
echo "raspd_root is now defined to be $raspd_root"
echo "path to conda is loaded"
echo "TRAPP installation is searched for at: $TRAPP"

