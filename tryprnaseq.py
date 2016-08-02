#!/usr/bin/env python3
#############################################################
#                       TrypRNAseq                          #
#          Automatic RNAseq raw file processing             #
#############################################################
#               created by Kevin Leiss                      #
#  individual parts by: Clementine Merce & Elisha Muchunga  #
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

import glob, os, datetime, itertools, subprocess, time, sys, importlib.util
from multiprocessing.dummy import Pool as ThreadPool
from threading import Thread

def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

script_path = getScriptPath()
sys.path.append(script_path)

import Rseq
#################################################################################

#from optparse import OptionParser

# Read in the variables, given by commandline options
options = Rseq.terminal_options(script_path)
if options.extension in ['fastq', 'fasta', 'gz']:
    ext = options.extension
    bow_indx = options.bow_index
    genome_gtf = options.gtf
    thread_no = options.threads
    exec_adapters = options.fastqc
    no_k = options.max_align

    exec_cutadapt = options.remove_adapters

    if(exec_cutadapt == 'y' and exec_adapters == 'n'):
        exec_adapters = 'y'
        print('\n\nCutadapt without FastQC not possible FastQC will be executed')


    site = options.adapter_site
    min_len = options.min_length
    adap_max = options.max_adapters

    files = glob.glob('./*.' + ext)
    # Generating a list of to be processed filenames:
    fnames = []
    for f in files:
        fname = f.split('.')[1]
        fname = fname[1:]
        fnames.append(fname)


    Rseq.print_line()

else:

    # Let the user specify the file-extension of the files to be processed
    ext = input('\nSpecify file extension of the raw-data (without \'.\', last extension [ie. txt.gz --> gz]): ')
    # Check whether the file extension is fastq or fasta.
    # This is important to prevent cutadapt from
    # prompting an error.
    while ext not in ['fastq', 'fasta', 'gz']:
        print('Please rename your reads to the right file-extension [gz (if compressed), fastq or fasta]')
        ext = input('Specify file extension: ')

    # Identifying the files with the previously given file extension:
    files = glob.glob('./*.' + ext)
    # Generating a list of to be processed filenames:
    fnames = []
    for f in files:
        fname = f.split('.')[1]
        fname = fname[1:]
        fnames.append(fname)


    Rseq.print_line()
    # Checking whether files are present with the specified extension
    while len(fnames) == 0:
        print('\n\nFiles can not be found, make sure you used the correct file-extension.')
        ext = input('File extension [gz, fasta, fastq]: ')
        files = glob.glob('./*.' + ext)
        fnames = []
        for f in files:
            fname = f.split('.')[1]
            fname = fname[1:]
            fnames.append(fname)


    # Ask if default parameters should be used
    print(('\n\nDefault-parameters: \n'
        '- Provided Bowtie Tryp. index \n'
        '- Provided .gtf file for read count \n'
        '- Remove all found adapters on both sides \n'
        '- Keep a minimal length of 30bp/read, discard all shorter \n'
        '- Number of threads = Will be asked for'))

    exec_default = input('Should the pipeline be executed with default parameters? y/n :')
    thread_no = input('How many numbers of threads should be used?: ')


    if exec_default in ['y', 'yes', 'Y']:
        bow_indx = script_path + '/bowtieindex/TbGenome'
        genome_gtf = script_path + '/Tb_cds.gtf'
        exec_adapters = 'y'
        exec_cutadapt = 'y'
        site = 'b'
        min_len = '30'
        adap_max = 'all'

# If not, let the user specify the parameters
    else:

        # Getting the bowtie2 index file
        bow_indx = input(
            'Please enter the path and the prefix (eg. Tbgenome) to your bowtie genome index: ')

        # Getting the gtf file for read counting
        genome_gtf = input('Path to genome gtf file: ')


        # Getting informations on how the data should be processed
        # If only the quality control of FastQC should be used,
        # but another list of adapters be removed, place a fasta file
        # with the adapters in the subfolder 'adapters'.
        # The file should have the same name as the file beeing processed with the extension:
        # _adapters.fasta
        exec_adapters = input(
            '\nShould a list of adapters be produced by FastQC? y/n: ')
        exec_cutadapt = input('\nShould the adapters be removed? y/n: ')

        if exec_cutadapt == 'y':
            if not os.path.exists('./adapters'):
                print('\n\nThe adapters-folder does not exist! Adapters will be generated...')
                exec_adapters = 'y'

            site = input(
                'Where are the adapters located? 3(a), 5(g) or both possible(b): ')
            min_len = input(
                'What is the minimal sequence length which should be kept?: ')
            adap_max = input(
                'How many adapters should be used for removal (the more the longer it takes)[int OR all]: ')


# Creating a log-file
settingslog = open('logfile.log', 'w')
settingslog.write('-----TrypRNAseq-----\n'
    'Pipeline started:\t' + str(datetime.datetime.now().date()) + '\t' + str(datetime.datetime.now().time()) + '\n'
    'Used pipeline settings:\n\n'
    'BowtieIndex:\t' + bow_indx + '\n'
    'GTF-File:\t' + genome_gtf + '\n'
    'Adapters removed?\t' + exec_cutadapt + '\n'
    'Minimal kept length:\t' + min_len + '\n'
    'Number of rem. adapters:\t' + adap_max)


Rseq.print_line()
# Unzipping gzipped files and adding correct extension

# This runs, if command line tool is used
# if options.extension in ['fastq', 'fasta', 'gz']:
#     if ext == 'gz':
#         os.system('mkdir gzipped_reads')
#         print('\nYour files are compressed. They will be decompressed.')
#         print('Extracting files')
#         ext = options.ext_unzip
#         pool = ThreadPool(int(thread_no))
#         pool.starmap(Rseq.gz_process, zip(files, itertools.repeat(ext)))
#         pool.close()
#         pool.join()
#
# # This starts if the interactive dialogue is used to set up the pipeline
# else:
if ext == 'gz':
    os.system('mkdir gzipped_reads')
    print('\nYour files are compressed. They will be decompressed.')
    if options.extension == 'gz': # checks whether commandline tool is used
        ext = options.ext_unzip
    else:  # or the interactive dialogue
        ext = input('Please specify the file extension of the decompressed file [fasta, fastq]: ')
    # Parallelized extraction and copying
    print('Extracting files')
    pool = ThreadPool(int(thread_no))
    pool.starmap(Rseq.gz_process, zip(files, itertools.repeat(ext)))
    pool.close()
    pool.join()

# Executing the FastQC algorithm
adap_set = adap_max
Rseq.print_line()
if exec_adapters in ['y', 'Y', 'yes']:
    # Analyzing the data with FastQC
    print('\nFastQC data analysis\n')
    Rseq.fqc([file_name + '.' + ext for file_name in fnames], thread_no)
    print('\nFastQC finished\n')

    # Generating the adapter list
    for fname in fnames:

        # Generating a adapters list
        # containing the 80 most abundant adapters
        print('\n\n' + fname + ' adapter list will be generated.')
        ex_adapters = Rseq.extract_adapters(fname)

        if adap_set == 'all':
            adap_max = len(ex_adapters)
        else:
            adap_max = int(adap_max)

        Rseq.make_ad_fasta(ex_adapters, fname, no_adapters=adap_max)
        print('\n\n' + fname + ' adapter list generated\n\n')


# Remove the adapters, stored in 'adapters'
Rseq.print_line()
if exec_cutadapt in ['y', 'Y', 'yes']:
    print('Cutadapt started')
    commands, summaries = Rseq.cutadapt(fnames, ext, site, min_len)
    processes = []
    for command in commands:
        processes.append(subprocess.Popen(command, shell=True))
    while 1:
        status = []
        status = [x.poll() for x in processes]
        if None not in status:
            break
        else:
            print(str(datetime.datetime.now().date()) + '\t' + str(datetime.datetime.now().time()) + ': Adapter removal running')
            time.sleep(5)

    processes = []
    print('Summarizing Cutadapt Results')
    for sum_cut in summaries:
        processes.append(subprocess.Popen(sum_cut, shell=True))
    while 1:
        status = []
        status = [x.poll() for x in processes]
        if None not in status:
            break
    print('Cutadapt done')

# Running bowtie either on the trimmed reads....
Rseq.print_line()
if exec_cutadapt in ['y', 'Y', 'yes']:
    for fname in fnames:

        fpath = './rm_adapt/' + fname + '/' + fname + '_processed.fastq'
        # Here bowtie is started using the processed data
        print('\n\nStarting alignment of ' + fname)
        Rseq.bowtie(fname, filepath=fpath, bow_index=bow_indx, no_threads = thread_no, k = no_k)
# ... or on the original reads
else:
    for fname in fnames:
        fname_ext = fname + '.' + ext
        fpath = './' + fname_ext

        # Here bowtie is started using the raw data, if no adapter removal was
        # done
        print('\n\nStarting alignment of ' + fname)
        Rseq.bowtie(filename=fname, filepath=fpath, bow_index=bow_indx, no_threads = thread_no, k = no_k)


for fname in fnames:
    fpath = './bowalign/' + fname + '_bow.sam'
    Rseq.print_line()
    print('\n\nPreparing the BAM files out of SAM files for \n' + fname)
    Rseq.sam_process(filename=fname, filepath=fpath, no_threads = thread_no)


for fname in fnames:
    fpath = './bam_files/' + fname + '_sorted.bam'
    Rseq.print_line()
    print('Generating index files for\n' + fname)
    Rseq.sam_index(fpath)

# Doing the read count
#for fname in fnames:
#    fpath = './bam_files/' + fname + '_sorted.bam'
#    Rseq.print_line()
#    print('Counting Reads for ' + fname)
#    Rseq.cds_only_counts(genome_gtf, fpath, fname)

fpaths = []
for fname in fnames:
    fpaths.append('./bam_files/' + fname + '_sorted.bam')
pool = ThreadPool(int(thread_no))
pool.starmap(Rseq.cds_only_counts, zip(itertools.repeat(genome_gtf), fpaths, fnames))
pool.close()
pool.join()

# Summarizing the generated read counts into one table
fpath = ''
reads_dicts = {}
for fname in fnames:
    fpath = './reads/' + fname + '_gene_read.txt'
    reads_dicts[fname] = Rseq.read_reads(fpath)

Rseq.create_reads_table(reads_dicts, 'reads')



settingslog.write('\n\nPipeline finished at:\t' + str(datetime.datetime.now().date()) + '\t' + str(datetime.datetime.now().time()))
print(3*'\n')
Rseq.print_line()
Rseq.print_line()
print('PIPELINE FINISHED')
Rseq.print_line()
Rseq.print_line()
print(3*'\n')
os.system('open ./reads/')
