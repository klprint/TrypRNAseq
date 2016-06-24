# TrypRNAseq
The TrypRNAseq pipeline is designed to take raw Illumina reads, does quality control, removal of overrepresented sequences and alignes the processed reads to a genome. In the end it counts the reads to a user-supplied .gtf file and produces a tab-sepparated file summarizing its results. Users should only specify in the beginning parameters which are asked for in a command-line dialogue.

## Information
- This code is still pre-alpha and is now under construction and cleanup!
- Till now, all reads will be aligned via bowtie2 up to 20-times (-k 20); in future releases this might be changed

## Dependencies
All dependencies need to be reachable via the command line.
- __Python version 3__
- [FastQC v0.11.5](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Cutadapt version 1.7.1](https://cutadapt.readthedocs.io/en/stable/)
- [bowtie2 version 2.0.0-beta7](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

## Usage
In the following, the usage of TrypRNAseq is described:


1. Clone the repository to your desired destination
2. Add the files containing the reads (can be gzipped, NO tarballs, make sure all have the same extension)
3. Start the pipeline using: python3 rnaseq_auto.py
4. Follow the instructions prompted in the terminal
5. If alignment against _Trypanosoma brucei_ TREU927 genome is intended, the included bowtie2 index can be used
6. If read counting should be done, using the coding sequences of the genes, the delivered GTF file (Tb_cds.gtf) can be used
7. After the pipeline finished, the folder with tab separated read counts will open automatically
8. For data analysis take care, that some genes are annotated by TrypDB with two, or more, CDS, therefore sum all reads up, originating in the same geneID

Now you can choose between default or user specified settings. Default settings are:

```
Default-parameters: 
   - Provided Bowtie Tryp. index 
   - Provided .gtf file for read count 
   - Remove all found adapters on both sides 
   - Keep a minimal length of 30bp/read, discard all shorter 
   - Number of threads = Will be asked for
```

## Licence
This pipeline was created by Kevin Leiss of the [Clayton lab](http://www.zmbh.uni-heidelberg.de/clayton/default.shtml) (ZMBH, Centre for Molecular Biology Heidelberg, Germany).


<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">TrypRNAseq</span> by <a xmlns:cc="http://creativecommons.org/ns#" href="https://github.com/klprint/TrypRNAseq" property="cc:attributionName" rel="cc:attributionURL">Kevin Leiss</a> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.