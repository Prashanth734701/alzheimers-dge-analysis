# Genome FASTA
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

# Annotation GTF
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz
gunzip gencode.v46.annotation.gtf.gz

# Build HISAT2 index
hisat2-build -p 8 GRCh38.primary_assembly.genome.fa GRCh38_index