#!/bin/bash
# Run in betaenv - Python v2.7
# ERROR! Python version 3.8 is not supported!
# Supported versions are 2.5, 2.6, 2.7, 3.3, 3.4, 3.5
# quast contig fix step

ProgDir=/home/akinya/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_lactucae/race_1/race_1_smartdenovo.dmo.lay.utg); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=assembly/SMARTdenovo/$Organism/$Strain/ncbi_edits
    sbatch $ProgDir/sub_quast.sh $Assembly $OutDir
  done


# Short cut to check contig number: less .fasta | grep '>' | wc -l
