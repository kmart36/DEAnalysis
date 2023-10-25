# DEAnalysis
### Differential Expression Analysis Pipeline for paired-end Illumina reads.

#### Programs used: 
FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ \
MultiQC: https://multiqc.info/ \
Trimmomatic: https://github.com/timflutre/trimmomatic \
Kraken2: https://ccb.jhu.edu/software/kraken2/ \
Trinity: https://github.com/trinityrnaseq/trinityrnaseq/wiki \
Salmon: https://combine-lab.github.io/salmon/ \
BUSCO: https://busco.ezlab.org/ \
cd-hit: https://sites.google.com/view/cd-hit/home?authuser=0 \

Important! Before running this pipeline, you need a metadata file for your raw data. This should be in .csv format and includes Locality, Tissue, and Sex for each sample. An example is included titled 'latra_metadata.csv'. 
