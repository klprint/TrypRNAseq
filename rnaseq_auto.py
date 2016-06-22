#############################################################
#                       RNAseq_auto                         #
#          Automatic RNAseq raw file processing             #
#############################################################
#               created by Kevin Leiss                      #
#          last updated on April 20th, 2016                 #
#############################################################

# Make sure you have installed:
# -fastqc
# -bowtie2
# -cutadapt
# -samtools

# Place the script in a folder with the raw data
# and make sure that only these data have the later specified
# file-extension.
# Then lean back and enjoy a coffee.. or three since 
# analysis will take a while.
#
# What you get in the end is a table consisting of 
# the gene-IDs and the associated reads.

import glob
import os
import Rseq


# Let the user specify the file-extension of the files to be processed

ext = raw_input('\nSpecify file extension of the raw-data (without .): ')

files = glob.glob('./*.' + ext)
bow_indx = raw_input(
    'Please enter the path and the prefix (eg. Tbgenome) to your bowtie genome index: ')

genome_gtf = raw_input('Path to genome gtf file: ')


# Getting informations on how the data should be processed
exec_adapters = raw_input(
    '\nShould a list of adapters be produced by FastQC? y/n: ')
exec_cutadat = raw_input('\nShould the adapters be removed? y/n: ')

if exec_cutadat == 'y':
    if not os.path.exists('./adapters'):
        print '\n\nThe adapters-folder does not exist! Adapters will be generated...'
        exec_adapters = 'y'

    site = raw_input(
        'Where are the adapters located? 3(a), 5(g) or both possible(b): ')
    min_len = raw_input(
        'What is the minimal sequence length which should be kept?: ')
    adap_max = raw_input(
        'How many adapters should be used for removal (the more the longer it takes)[int OR all]: ')



# Generating a list of to be processed files:
fnames = []
fname_ext = ''
for f in files:
    fname = f.split('.')[1]
    fname = fname[1:]
    fname_ext = fname + '.' + ext
    fnames.append(fname)

if exec_adapters == 'y':
    for fname in fnames:

        # Analyzing the data with FastQC
        print '\nFastQC data analysis\n'
        Rseq.fqc(fname_ext)
        print '\nFastQC finished\n'

        # Generating a adapters list
        # containing the 80 most abundant adapters
        print '\n\n' + fname + ' adapter list will be generated.'
        ex_adapters = Rseq.extract_adapters(fname)

        if adap_max == 'all':
            adap_max = len(ex_adapters)
        else:
            adap_max = int(adap_max)

        Rseq.make_ad_fasta(ex_adapters, fname, no_adapters=adap_max)
        print '\n\n' + fname + ' adapter list generated\n\n'


if exec_cutadat == 'y':
    for fname in fnames:
        # This function removes the adaptors found before
        Rseq.cutadapt(fname, ext, site, seq_min_len=min_len)


if exec_cutadat == 'y':
    for fname in fnames:

        fpath = './rm_adapt/' + fname + '/' + fname + '_processed.fastq'
        # Here bowtie is started using the processed data
        Rseq.bowtie(fname, filepath=fpath, bow_index=bow_indx)
else:
    for fname in fnames:
        fpath = './' + fname + '.' + fname_ext

        # Here bowtie is started using the raw data, if no adapter removal was
        # done
        Rseq.bowtie(filename=fname, filepath=fpath, bow_index=bow_indx)


for fname in fnames:
    fpath = './bowalign/' + fname + '_bow.sam'

    print '\n\nPreparing the BAM files out of SAM files for \n' + fname
    Rseq.sam_process(filename=fname, filepath=fpath)


for fname in fnames:
    fpath = './bam_files/' + fname + '_sorted.bam'

    print 'Generating index files for\n' + fname
    Rseq.sam_index(fpath)


for fname in fnames:
    fpath = './bam_files/' + fname + '_sorted.bam'

    Rseq.cds_only_counts(genome_gtf, fpath, fname)


