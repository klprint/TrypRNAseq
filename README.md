# TrypRNAseq
The TrypRNAseq pipeline is designed to take raw Illumina reads, does quality control, removal of overrepresented sequences and alignes the processed reads to a genome. In the end it counts the reads to a user-supplied .gtf file and produces a tab-sepparated file summarizing its results. Users should only specify in the beginning parameters which are asked for in a command-line dialogue.

## Dependencies
All dependencies need to be reachable via the command line.
- [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

## Information
This code is still pre-alpha and is now under construction and cleanup!