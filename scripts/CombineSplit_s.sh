#Output files: DSA1*_interpro.gff3 , DSA*_interproscan.tsv , DSA1*_interproscan.xml
#outdir gene_pred/interproscan/F.oxysporum_fsp_fragariae/DSA*

ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  for Proteins in $(ls gene_pred/codingquary/F.oxysporum_fsp_fragariae/DSA15_041/final/final_genes_appended_renamed.pep.fasta); do
    Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    echo $Strain
    InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
    $ProgDir/append_interpro.sh $Proteins $InterProRaw
  done
