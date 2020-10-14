# Convert Salmon quasi-quanitifcations to gene counts using an awk script:


mkdir -p RNAseq_analysis/salmon/F.oxysporum_fsp_fragariae/DSA14_003/DeSeq2
#This command creates a two column file with transcript_id and gene_id.
for File in $(ls RNAseq_analysis/salmon/F.oxysporum_fsp_fragariae/DSA14_003/*/quant.sf | head -n1); do
  cat $File | awk -F"\t" '{c=$1;sub(".t.*","",$1);print c,$1}' OFS="\t" > RNAseq_analysis/salmon/F.oxysporum_fsp_fragariae/DSA14_003/DeSeq2/trans2gene.txt
done

#This command creates a two column file with transcript_id.
#for File in $(ls alignment/salmon/*/Hg199/*/*/*/quant.sf | head -n1); do
#cat $File | awk -F"\t" '{c=$1;sub("\*","",$1);print c,$1}' OFS="\t" > alignment/salmon/N.ditissima/Hg199/DeSeq2/trans2gene3.txt
#done

# Put files in a convenient location for DeSeq.
# Analysis was not performed on Apple control samples.

for File in $(ls RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/*/quant.sf); do
  Prefix=$(echo $File | cut -f5 -d '/' --output-delimiter '_')
  mkdir -p RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/$Prefix
  cp $PWD/$File RNAseq_analysis/salmon_vAG/F.venenatum/WT_minion/DeSeq2/$Prefix/quant.sf
  # rm alignment/salmon/DeSeq2/$Prefix/quant.sf
done
