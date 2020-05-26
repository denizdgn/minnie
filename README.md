# minnie
A structural ensemble analysis package for extracting molecular interaction fingerprints

minnie functions under the umbrella of anaconda3 (with Python 3.6). So, before everything you will need to have anaconda3 installed.

## Clone the repository
```
git clone https://github.com/CSB-KaracaLab/minnie.git
```
OR
```
wget https://github.com/CSB-KaracaLab/minnie/archive/master.zip
```

## Set up the right environment for minnie

```
cd minnie

bash setup.sh

```
setup.sh script installs the required dependencies together with the first version of interfacea (https://zenodo.org/badge/latestdoi/136096537). 
It finally calls execute.sh script, which sets proper aliases and updates your .bashrc file. 

## Activate minnie

When the environment is set
```
source ~/.bashrc
```

## Run minnie 

These commands should be run in your project folder where you ensembles are located. 

To split trajectories into single frames where your project ids are sox4 & sox18
```
minnie splitpdbs -cn sox4 sox18 -p sox4.pdb sox18.pdb
```

This part calculates all types of inter-monomer interactions
```
minnie findbonds -cn sox4 -p sox4/02_frames/* -i all

minnie findbonds -cn sox18 -p sox18/02_frames/* -i all
```

Apply critical interaction filter
```
minnie timefilter -f sox4/03_interfacea_results/*/sox4_merged_*.csv -cn sox4 --per 25

minnie timefilter -f sox18/03_interfacea_results/*/sox18_merged_*.csv -cn sox18 --per 25
```

Calculate common and distinct interactions in two cases
```
minnie compareCX -cn sox4 sox18 --per 25
```

aaaaand build graphs!
```
minnie graph -cn 'sox4' 'sox18' --per 25 -i all -b protein-dna -s specific

minnie graph -cn 'sox4' 'sox18' --per 25 -i all -b protein-dna -s common
```
All of these commands are also provided in pipeline.sh

## Troubleshoot
If you would need to remove minnie from your conda setup
```
conda env remove -n  minnie
```
Also remove the minnie-related lines from your .bashrc

## Contact
ezgi.karaca@ibg.edu.tr
