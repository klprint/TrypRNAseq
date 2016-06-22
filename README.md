# TrypRNAseq
The TrypRNAseq pipeline is designed to take raw Illumina reads, does quality control, removal of overrepresented sequences and alignes the processed reads to a genome. In the end it counts the reads to a user-supplied .gtf file and produces a tab-sepparated file summarizing its results. Users should only specify in the beginning parameters which are asked for in a command-line dialogue.

## Information
This code is still pre-alpha and is now under construction and cleanup!

## Dependencies
All dependencies need to be reachable via the command line.
- [FastQC v0.11.5](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Cutadapt version 1.7.1](https://cutadapt.readthedocs.io/en/stable/)
- [bowtie2 version 2.0.0-beta7](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)


## Licence
This pipeline was created by Kevin Leiss of the [Clayton lab](http://www.zmbh.uni-heidelberg.de/clayton/default.shtml) (ZMBH, Centre for Molecular Biology Heidelberg, Germany).
<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">TrypRNAseq</span> by <a xmlns:cc="http://creativecommons.org/ns#" href="https://github.com/klprint/TrypRNAseq" property="cc:attributionName" rel="cc:attributionURL">Kevin Leiss</a> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.