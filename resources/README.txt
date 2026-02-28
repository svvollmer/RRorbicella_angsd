Reference genome and annotation files live in:
  /home/minipc/coral-angsd-pipeline/reference/

These files are too large for git and are excluded from version control.

Expected contents:
  reference/GCF_021335395.2_Apul_2.0_genomic.fna          - A. palmata genome (FASTA)
  reference/GCF_021335395.2_Apul_2.0_genomic.fna.fai      - samtools faidx index
  reference/GCF_021335395.2_Apul_2.0_genomic.fna.amb      - BWA index files
  reference/GCF_021335395.2_Apul_2.0_genomic.fna.ann
  reference/GCF_021335395.2_Apul_2.0_genomic.fna.bwt
  reference/GCF_021335395.2_Apul_2.0_genomic.fna.pac
  reference/GCF_021335395.2_Apul_2.0_genomic.fna.sa
  reference/GCF_021335395.2_Apul_2.0_genomic.gff          - Gene annotation (optional)

To download the reference genome:
  datasets download genome accession GCF_021335395.2 --include genome,gff3
