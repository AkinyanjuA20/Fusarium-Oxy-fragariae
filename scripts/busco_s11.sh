#Edited paths from Andy's direc to mine
#Updated entire script using https://github.com/harrisonlab/bioinformatics_tools/blob/master/Gene_prediction/README.md
#Does not have quast segment in script like Andy's pipeline
#Look into BuscoDB direc - directory exists
#Run in conda env - Repenv

for Assembly in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_unmasked.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
ProgDir=/home/akinya/git_repos/tools/seq_tools/busco
BuscoDB=$(ls -d /projects/oldhome/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
sbatch $ProgDir/sub_busco.sh $Assembly $BuscoDB $OutDir
done

