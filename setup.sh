#!/usr/bin/env bash

# Set the right environment for minnie

conda env create -f environment.yml
conda activate minnie
cd interfacea
python setup.py build
python setup.py install
cd ..
bash execute.sh