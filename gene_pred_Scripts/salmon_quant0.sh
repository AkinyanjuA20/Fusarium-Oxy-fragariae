# Run in Salmon env
# This script simply loops through each sample and invokes salmon
# -i - where to find the index
# -l A - automatically determines the library type of the sequencing reads
# -1 and -2 - where to find the left and right reads for this sample
# -p 8 - threads count
# -o output dir

#!/bin/bash
for fn in data/DRR0161{25..40};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i athal_index -l A \
         -1 ${fn}/${samp}_1.fastq.gz \
         -2 ${fn}/${samp}_2.fastq.gz \
         -p 8 --validateMappings -o quants/${samp}_quant
done
