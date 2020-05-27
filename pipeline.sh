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


#To visualize minnie results with Pymol!
for projectID in  sox4 sox18
do
  cd "$projectID/05_compare_complex/25_freq_filtered/25_freq_perres/complex_specific"
  sed 's/,/ /g' "$projectID_hbonds_compared_25_perres.csv" | grep "protein-dna" | grep ":" | awk '{printf "show sticks, (resi %s and chain %s) + (resi %s and chain %s) \n", $6,$2,$7,$3}' | sort -u > pymol_hbonds.pml
  sed 's/,/ /g' "$projectID_hbonds_compared_25_perres.csv" | grep "protein-dna" | grep ":" | awk '{printf "distance i. %s and n. %s and chain %s, i. %s and n. %s and chain %s\n", $6,$8,$2,$7,$9,$3}' | sort -u | sed 's/OP1/O1P/g' | sed 's/OP2/O2P/g' >> pymol_hbonds.pml


  sed 's/,/ /g' "$projectID_hydrophobic_compared_25_perres.csv" | grep "protein-dna" | grep ":" | awk '{printf "show sticks, (resi %s and chain %s) + (resi %s and chain %s) \n", $6,$2,$7,$3}' | sort -u > pymol_hydrophobic.pml
  sed 's/,/ /g' "$projectID_hydrophobic_compared_25_perres.csv" | grep "protein-dna" | grep ":" | awk '{printf "distance i. %s and n. %s and chain %s, i. %s and n. %s and chain %s\n", $6,$8,$2,$7,$9,$3}' | sort -u | sed 's/OP1/O1P/g' | sed 's/OP2/O2P/g' >> pymol_hydrophobic.pml


  sed 's/,/ /g' "$projectID_ionic_compared_25_perres.csv" | grep "protein-dna" | grep ":" | awk '{printf "show sticks, (resi %s and chain %s) + (resi %s and chain %s) \n", $6,$2,$7,$3}' | sort -u > pymol_ionic.pml
  sed 's/,/ /g' "$projectID_ionic_compared_25_perres.csv" | grep "protein-dna" | grep ":" | awk '{printf "distance i. %s and n. %s and chain %s, i. %s and n. %s and chain %s\n", $6,$8,$2,$7,$9,$3}' | sort -u | sed 's/OP1/O1P/g' | sed 's/OP2/O2P/g' >> pymol_ionic.pml

  sed 's/,/ /g' "$projectID_ring_stacking_compared_25_perres.csv" | grep "protein-dna" | grep ":" | awk '{printf "show sticks, (resi %s and chain %s) + (resi %s and chain %s) \n", $6,$2,$7,$3}' | sort -u > pymol_ring_stacking.pml
  sed 's/,/ /g' "$projectID_ring_stacking_compared_25_perres.csv" | grep "protein-dna" | grep ":" | awk '{printf "distance i. %s and n. %s and chain %s, i. %s and n. %s and chain %s\n", $6,$8,$2,$7,$9,$3}' | sort -u | sed 's/OP1/O1P/g' | sed 's/OP2/O2P/g' >> pymol_ring_stacking.pml
  cd -
done
