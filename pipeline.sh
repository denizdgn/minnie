#!/usr/bin/env bash

#Splitting trajectories into single frames
minnie splitpdbs -cn sox4 sox18 -p sox4.pdb sox18.pdb

# This part calculates all types of inter-monomer interactions
minnie findbonds -cn sox4 -p sox4/02_frames/* -i all
minnie findbonds -cn sox18 -p sox18/02_frames/* -i all

# Apply critical interaction filter
minnie timefilter -f sox4/03_interfacea_results/*/sox4_merged_*.csv -cn sox4 --per 25
minnie timefilter -f sox18/03_interfacea_results/*/sox18_merged_*.csv -cn sox18 --per 25

# Calculate common and distinct interactions in two cases
minnie compareCX -cn sox4 sox18 --per 25

# aaaaand graphs!
minnie graph -cn 'sox4' 'sox18' --per 25 -i all -b protein-dna -s specific
minnie graph -cn 'sox4' 'sox18' --per 25 -i all -b protein-dna -s common