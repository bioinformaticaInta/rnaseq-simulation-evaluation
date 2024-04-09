#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import commands
from popen2 import popen2
import time
import numpy
import argparse
import random

################Input parsing arguments and parameters##########################################

parser = argparse.ArgumentParser(description='Script to simulate N samples of Rna-seq reads using art_illumina software.')
parser.add_argument("-gn", "--genes_number", dest='gene_N', metavar='<int>', default=10000, help="Number of genes to simulate reads (Default 10000 or total amount of genes in the organism).")
parser.add_argument("-sp", "--max_variants", dest='iso_max', metavar='<int>', default=1, help="Maximum number of splice variant for each gene (Default 1).")
parser.add_argument("-moi", "--more_one_isoform", dest='more_one', help="Only use genes with 2 or more isoforms.", action='store_true')
parser.add_argument("-es", "--error_shift", dest='errorsh', metavar='<int>', default=0, help="Value for the error shift. A value of 1000 is aproximately zero sequencing erros (Default=0)")
parser.add_argument("-sr", "--simulate_reads", dest='sim_reads', metavar='<dir>', default=None, help="Only simulate reads, the parameter need a directory with one fasta file for each transcript. (Default None)")
parser.add_argument("-em", "--expression_mean", dest='mean', metavar='<int>', default=4, help="Mean value for lognormal distribution (Default 4).")
parser.add_argument("-estd", "--expression_std", dest='std', metavar='<int>', default=1, help="Standard deviation for lognormal distribution (Default 1).")
parser.add_argument("-se", "--single_end", help="Simulate single end reads (Default is pair end with fragment of 400 bps).", action='store_true')
parser.add_argument("-pe", "--pair_end", dest='frag_len', metavar='<int>', default=400, help="Pair end fragment length (Default 400).")
parser.add_argument("-l", "--read_length", dest='read_len', metavar='<int>', default=150, help="Length of the reads (Default 150).")
parser.add_argument("-fna", "--fna_files", metavar='<file1>,<file2>', help="List of fna files for the reference transcripts.")
parser.add_argument("-gtf", "--gtf_file", metavar='<file>', help="File for reference transcripts descriptions in gtf format.")
parser.add_argument("-dup", "--dup_file", metavar='<file>', help="File for duplicate genes in reference (list).")
parser.add_argument("-d", "--output_directory", dest='out_dir', metavar='<dir>', default='rnaseq_simulation_files', help="Output directory with fastq files as results. (Default ./rnaseq_simulation_files).")
parser.add_argument("-n", "--fastq_files_number", dest='fq_N', metavar='<int>', default=10, help="Amount of generated fastq files (Default 10).")
parser.add_argument("-p", "--prefix", metavar='<prefix>', default='reads', help="Prefix for the output files, fastq and expression lists. Final names: <prefix>_N_R1.fastq <prefix>_N_R2.fastq and <prefix>_N.exp (Default reads).")

args = parser.parse_args()


################################################################################################

def waiting_empty_queue(message):
    cant_queue = commands.getoutput('qstat | wc -l')
    while int(cant_queue) > 0:
        print '%s: %s' % (message, cant_queue)
        time.sleep(30) #Waiting 30 seconds
        cant_queue = commands.getoutput('qstat | wc -l')
    return
    
if (args.fna_files and args.gtf_file) or args.sim_reads:
    #Create output directory or empty if exists
    os.system('if [ -d %s ]; then rm -r %s/*; else mkdir %s; fi' % (args.out_dir, args.out_dir, args.out_dir))
    #Write information in log file
    logfile = '%s/simulation.log' % (args.out_dir)
    log_file = open(logfile, 'w')
    log_file.write('Gene number: %s\nMax isoforms number: %s\nExpression mean: %s\nExpression std: %s\n' % (args.gene_N, args.iso_max, args.mean, args.std))
    log_file.write('Only use genes with 2  or more isoforms: %s\n' % ('True' if args.more_one else 'False'))
    log_file.write('Deduplicated genes: %s\n' % ('True' if args.dup_file else 'False'))
    log_file.write('Reads type: %s\n' % ('Single end' if args.single_end else 'Pair end'))
    log_file.write('Read length: %s\nFna files: %s\nGtf file: %s\nFastq files number: %s\n' % (args.read_len, args.fna_files, args.gtf_file, args.fq_N))
    log_file.close()

    if not args.sim_reads:
        #Merge all fna files in one file
        for fna_file in args.fna_files.split(','):
            os.system('cat %s >> %s/%s_complete_reference.fna' % (fna_file, args.out_dir, args.prefix))
       
        # Processing gtf file to create a dictionary with all transcripts for each gene 
        all_gene_dic = {}
        gtf_file = open(args.gtf_file)
        for line in gtf_file:
            line = line[:-1].split('\t')
            if '#' in line[0]: continue
            if line[2] == "exon":
                data = line[8].split(';')
                gene_name = transcript_name = ''
                for field in data:
                    fieldv = field.split('"')
                    if "gene_id" in fieldv[0]:
                        gene_name = fieldv[1]
                    if "transcript_id" in fieldv[0]:
                        transcript_name = fieldv[1]
                try:
                    all_gene_dic[gene_name].keys()
                except:
                    all_gene_dic[gene_name] = {}
                try:
                    all_gene_dic[gene_name][transcript_name] = all_gene_dic[gene_name][transcript_name] + int(line[4]) - int(line[3]) + 1
                except:
                    all_gene_dic[gene_name][transcript_name] = int(line[4]) - int(line[3]) + 1
        gtf_file.close()
	
        # Filtering the transcripts with length less than read lenght
        aux_dic = all_gene_dic
        genes = aux_dic.keys()
        if args.single_end:
            len_limit = int(args.read_len)
        else:
            len_limit = int(args.frag_len)
            
        for gene in genes:
            for trans in aux_dic[gene].keys():
                if aux_dic[gene][trans] < len_limit:
                    del all_gene_dic[gene][trans]
            if len(all_gene_dic[gene]) == 0:
                del all_gene_dic[gene]
        del aux_dic
        
        # Processing duplicated genes file and only genes with 2 or more isoforms filter
        # Filter duplicated genes
        gene_dic = {}
        if args.dup_file:
            dedup_gene_dic = {}
            dup_list = []
            dup_file = open(args.dup_file)
            for line in dup_file:
                line = line[:-1]
                dup_list.append(line)
            dup_file.close()
            for gene in all_gene_dic:
	        if gene not in dup_list:
		    dedup_gene_dic[gene] = all_gene_dic[gene]
 	    # Filter genes with less 1 isoform
	    if args.more_one:
                for gene in dedup_gene_dic:
	            if len(dedup_gene_dic[gene].keys()) > 1:
	                gene_dic[gene] = dedup_gene_dic[gene]
            else:
                gene_dic = dedup_gene_dic
        # Filter genes with less 1 isoform 
        elif args.more_one:
            for gene in all_gene_dic:
	        if len(all_gene_dic[gene].keys()) > 1:
	            gene_dic[gene] = all_gene_dic[gene]
	else:	
            gene_dic = all_gene_dic
        
        print len(gene_dic) 
        #Formating file to extract sequences after
        os.system('makeblastdb -in %s/%s_complete_reference.fna -dbtype nucl -parse_seqids' % (args.out_dir, args.prefix))
        
        #Generate separate fna files for each transcript
        for sample_N in range(1,int(args.fq_N)+1):
            os.system('mkdir %s/transcripts_fna_%s' % (args.out_dir, sample_N))
            transcript_list = []
            work_dir = os.getcwd()
            counter = 0
            transfile = '%s/gene_transcript_%s.list' % (args.out_dir, sample_N)
            trans_file = open(transfile, 'w')
            if int(len(gene_dic)) < int(args.gene_N): args.gene_N = len(gene_dic)
            for gene in random.sample(gene_dic.keys(), int(args.gene_N)):
                transcript_sorted = sorted(gene_dic[gene], key=gene_dic[gene].get, reverse=True)
                for transcript in transcript_sorted[0:int(args.iso_max)]:
                    if counter == 500:
                        counter = 0
                        waiting_empty_queue('Sequence to extract')
                    counter += 1
                    trans_file.write('%s,%s\n' % (gene, transcript))
                    transcript_list.append(transcript)
                    output, input = popen2('qsub')
                    input.write("cd %s; blastdbcmd -db %s/%s_complete_reference.fna -dbtype nucl -entry %s -out %s/transcripts_fna_%s/%s.fasta; sed -i 's/lcl|//' %s/transcripts_fna_%s/%s.fasta" % (work_dir, args.out_dir, args.prefix, transcript, args.out_dir, sample_N, transcript, args.out_dir, sample_N, transcript))
                    input.close()
            waiting_empty_queue('Sequence to extract')
            trans_file.close()    
    
            #Create files with expression levels and reads (fastq)
            work_dir = os.getcwd()
            counter = 0
            expfile = '%s/%s_%s.exp' % (args.out_dir, args.prefix, sample_N)
            exp_file = open(expfile, 'w')
            for transcript in transcript_list:
                if counter == 500:
                    counter = 0
                    waiting_empty_queue('Sequences to simulate reads')
                counter += 1
                expression_value = numpy.random.lognormal(mean=args.mean,sigma=args.std,size=1)[0]
                #Write expression values in list file
                exp_file.write('%s,%s\n' % (transcript,expression_value))
                #Simulate reads
                if args.single_end:
                    art_command = 'cd %s; art_illumina -na -qs %s -i %s/transcripts_fna_%s/%s.fasta -l %s -f %s -o %s' % (work_dir, args.errorsh, args.out_dir, sample_N, transcript, args.read_len, expression_value, transcript)
                else:
                    art_command = 'cd %s; art_illumina -na -qs %s -qs2 %s -i %s/transcripts_fna_%s/%s.fasta -p -l %s -f %s -m %s -s 200 -o %s_' % (work_dir, args.errorsh, args.errorsh, args.out_dir, sample_N, transcript, args.read_len, expression_value, args.frag_len, transcript)
            
                output, input = popen2('qsub')
                input.write(art_command)
                input.close() 
        
            exp_file.close()
            waiting_empty_queue('Sequences to simulate reads')

            #Merge fastq files
            if args.single_end:
                os.system('cat *.fq > %s/%s_%s.fastq' % (args.out_dir, args.prefix, sample_N))
                os.system('rm *.fq')
            else:
                os.system('cat *_1.fq > %s/%s_%s_1.fastq' % (args.out_dir, args.prefix, sample_N))
                os.system('cat *_2.fq > %s/%s_%s_2.fastq' % (args.out_dir, args.prefix, sample_N))
                os.system('rm *.fq')
            # Delete all individual fasta files for transcripts and create the reference file
            os.system('cat %s/transcripts_fna_%s/* > %s/all_reference_transcripts_%s.fna' % (args.out_dir, sample_N, args.out_dir, sample_N))
            os.system('rm -r %s/transcripts_fna_%s' % (args.out_dir, sample_N))

        # Delete file with complete reference
        os.system('rm %s/%s_complete_reference.fna*' % (args.out_dir, args.prefix))
