Script for importing FASTQ files consisting of multiple sequence sets with 4 lines:
- Line 1: Sample name/number/description
- Line 2: Raw sequence
- Line 3: Empty (gets skipped in this script)
- Line 4: Quality scores (encoded)

Features:
- Import of files in the above described format
- Option to output total "G" and "C" content in each data set
- Option to output content of defined single base or sequence in each data set
- Option to output data as FASTA formated file
- Output of average quality score for each data set
- Output of sample description for each data set
- Output of triplet converted to amino acid sequence for first sequence in file