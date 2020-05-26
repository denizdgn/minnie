#!/usr/bin/env bash

shopt -s expand_aliases
path_to_scripts=$(pwd)
chmod +x minnie.py


source ~/anaconda3/etc/profile.d/conda.sh

echo '##### Required for minnie #####' >> ~/.bashrc
echo "conda activate minnie" >> ~/.bashrc
echo "alias '"'minnie=python "'$path_to_scripts/minnie.py'"  '"' " >> ~/.bashrc
echo '##### Required for minnie #####' >> ~/.bashrc

source ~/.bashrc