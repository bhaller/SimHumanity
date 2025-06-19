This "mtDNA info" folder contains information about the genetic structure of the mitochondrial "chromosome".

This data was obtained from the NIH National Library of Medicine, NCBI Reference Sequence: NC_012920.1,
which is currently located at https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1

The total length, 16569 bp, used in the model is also taken from there.

This folder contains:

- this file
- chrMT_genomic_elements.txt - the coding sequence for the mtDNA
- mtDNA.FASTA : the original FASTA file for NC_012920.1
- mtDNA_MOD.FASTA : the modified FASTA, changing one N nucleotide to A

The mtDNA_MOD.FASTA file is a modification of mtDNA.FASTA that changes one nucleotide given as 'N' to 'A' so that SLiM is happy with it.

The chrMT_genomic_elements.txt file has columns (getype-id) (start) (end) in zero-based positions.  CDS segments are coded as g1 (with a 1), all other ranges in the chromosome are coded as g0 (with a 0).  Some CDS segments abutted or overlapped; these have simply been merged together.  One CDS segment is on the complementary strand; its range on the main strand is used, since SLiM does not model the two strands separately.
