#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import commands
import time
import numpy as np
import argparse
import random
import matplotlib.pyplot as plt

################Input parsing arguments and parameters##########################################

parser = argparse.ArgumentParser(description='Script to calculate statistics over the contigs result of RNA-Seq assembly based on comparation with the reference (for example quimera identification) and also without reference (for example N50).')
parser.add_argument("-as", "--assembly_software", dest='a_soft', metavar='<software_name>',choices=['oases', 'trinity', 'abyss', 'soap','idba','bridger'], default=None, help="Assembly software using in the asembly proccess. Options available: oases, trinity, abyss, soap, idba, bridger (Default None).")
parser.add_argument("-a", "--assembly_file", metavar='<file>', required=True, help="Fasta file with contigs of assembly results.")
parser.add_argument("-r", "--read_file", metavar='<file>', required=True, help="Fasq file with raw read.")
parser.add_argument("-ref", "--reference_file", dest='ref_file', metavar='<file>', required=True, help="Fasta file with reference transcripts.")
parser.add_argument("-t", "--gene_transcript", dest='gene_trans_file', metavar='<file>', required=True, help="File contains relationships between gene and transcripts on reference in csv format.")
parser.add_argument("-e", "--expression_file", dest='exp_file', metavar='<file>', required=True, help="File contains the expression value of each reference transcipt in csv format.")
parser.add_argument("-p", "--prefix_name", dest='prefix', metavar='<file>', default="output", help="Prefix for the output files with contig names, length, gc content, etc.")
parser.add_argument("-l", "--read_length", dest='read_len', metavar='<int>', default=150, help="Read length.")

args = parser.parse_args()

################################################################################################

def calculate_Nxx(list_length, xx):
    partial_sum = 0
    xx_length = sum(list_length)*(float(xx)/float(100))
    for length in list_length:
        partial_sum += length
        if partial_sum >= xx_length:
            return length

def load_in_ass_gapped_dic(dic, tname, cov, covi, gap_p, hname, bsize, bpos):
    if tname in dic:
#        if cov < float(1):
#             if cov > dic[tname][0]:
#                 dic[tname] = [cov, covi, gap_p, hname, bsize, bpos]
#        else:
         if covi > dic[tname][1]:
             dic[tname] = [cov, covi, gap_p, hname, bsize, bpos]
    else:
        dic[tname] = [cov, covi, gap_p, hname, bsize, bpos]

################################################################################################

if args.assembly_file:
    
    #######################################################################################################################
    #STATITICS VALUES INDEPENDENT OF THE REFERENCE
    #MEAN AND MEDIAN CONTIG LENGTH
    #N20, N50 AND N80
    #######################################################################################################################

    #Open output file
    out_file = open('%s_all_stats.txt' % (args.prefix), 'w')
    
    os.system('infoseq -nodescription -nousa -nodatabase -noaccession -nogi -noseqversion -notype -noorganism -noheading -sequence %s -outfile %s.infoseq' % (args.assembly_file, args.prefix))
    info_file = open('%s.infoseq' % args.prefix, 'r')
    info_file_list = info_file.readlines()
    info_file.close()
    contig_lengths = map(lambda x: int(x.split()[1]), info_file_list)
    contig_gcs = map(lambda x: float(x.split()[2]), info_file_list)
    mean_length = np.mean(contig_lengths)
    median_length = np.median(contig_lengths)
    contig_lengths.sort(reverse=True)
    N20 = calculate_Nxx(contig_lengths, 20)
    N50 = calculate_Nxx(contig_lengths, 50)
    N80 = calculate_Nxx(contig_lengths, 80)

    out_file.write("Mean of contig length = %s\n" % mean_length)
    out_file.write("Median of contig length = %s\n" % median_length)
    out_file.write("N20 = %s\n" % N20)
    out_file.write("N50 = %s\n" % N50)
    out_file.write("N80 = %s\n" % N80)

    #######################################################################################################################
    #BLAT alignments
    #######################################################################################################################
    
    os.system("blat %s /dev/null /dev/null -tileSize=11 -makeOoc=%s.ooc -repMatch=1024" % (args.assembly_file, args.assembly_file))
    os.system("blat %s /dev/null /dev/null -tileSize=11 -makeOoc=%s.ooc -repMatch=1024" % (args.ref_file, args.ref_file))
    os.system("blat -minIdentity=95 %s %s -ooc=%s.ooc %s_reference_to_assembly.psl" % (args.assembly_file, args.ref_file, args.assembly_file, args.prefix))
    os.system("blat -minIdentity=95 %s %s -ooc=%s.ooc %s_assembly_to_reference.psl" % (args.ref_file, args.assembly_file, args.ref_file, args.prefix))
    os.system("blat -minIdentity=95 %s %s -ooc=%s.ooc %s_assembly_to_assembly.psl" % (args.assembly_file, args.assembly_file, args.assembly_file, args.prefix))
     
    #######################################################################################################################
    #REF-EVAL STATISTICS
    #######################################################################################################################
    
#    read_number = commands.getoutput("wc -l %s | cut -d ' ' -f 1" % (args.read_file))
#    read_number = int(read_number) / 4 
#    os.system("rsem-prepare-reference --no-polyA %s assembly_ref" % (args.assembly_file))
#    os.system("rsem-prepare-reference --no-polyA %s reference_ref" % (args.ref_file))
#    os.system("rsem-calculate-expression --no-bam-output %s assembly_ref assembly_expr" % (args.read_file))
#    os.system("rsem-calculate-expression --no-bam-output %s reference_ref reference_expr" % (args.read_file))
#    os.system("ref-eval --scores=nucl,contig,kmer,kc --weighted=both --A-seqs %s --B-seqs %s --A-to-B %s_assembly_to_reference.psl --B-to-A %s_reference_to_assembly.psl --A-expr assembly_expr.isoforms.results --B-expr reference_expr.isoforms.results --readlen %s --num-reads %s --kmerlen 31 > %s_ref_eval_scores.txt" % (args.assembly_file, args.ref_file, args.prefix, args.prefix, args.read_len, read_number, args.prefix))

    #######################################################################################################################
    #Calculate amounts of reference and assembled transcripts 
    #######################################################################################################################
    
    all_ass_trans = commands.getoutput("fgrep '>' %s | cut -d ' ' -f 1 | sed 's/>//g'" % (args.assembly_file)).split('\n')
    all_ref_trans = commands.getoutput("fgrep '>' %s | cut -d ' ' -f 1 | sed 's/>//g'" % (args.ref_file)).split('\n')
    ass_trans_N = float(len(all_ass_trans))
    ref_trans_N = float(len(all_ref_trans))
    ass_reftrans_relation = ass_trans_N / ref_trans_N    

    out_file.write("Number of reference transcripts = %s\n" % ref_trans_N)
    out_file.write("Number of assembly transcripts = %s\n" % ass_trans_N)
    out_file.write("Assembly transcripts / Reference transcripts = %s\n" % ass_reftrans_relation)
    
    #######################################################################################################################
    #FULL LENGTH RECONSTRUCTED TRANSCRIPTS (99, 95, 90, 80 AND 70%)
    #FALSE POSITIVE RATE
    #NUCLEOTIDE SENSITIVITY AND ESPECIFICITY
    #######################################################################################################################

    ass_to_ref = open('%s_assembly_to_reference.psl' % args.prefix, 'r')
    
    full_rec99_dic = {}
    full_rec95_dic = {}
    full_rec90_dic = {}
    full_rec80_dic = {}
    full_rec70_dic = {}
    ass_gapped_dic = {}
    false_pos_list = []
    true_pos_list = []
    trans_with_hit_list = []
    true_max_match_len_dic = {}
    true_match_dic = {}

    line_flag=0
    #Parsing assembly to reference psl file
    for line in ass_to_ref:
        line_flag+=1
        #No work with first 5 lines of header file
        if line_flag>5:
            line = line.split()

            #Align length adding matches, missmatches, Ns  and indels
            indel_N = max(float(line[5]), float(line[7]))
            align_len_without_indels = float(line[0]) + float(line[1]) + float(line[2])
            align_len = align_len_without_indels + indel_N
            indel_p = indel_N / align_len
            
            #Calculate reference and assembly coverage and indels percentage
            ref_cov = align_len / float(line[14])
            ass_cov = align_len / float(line[10])
            ass_cov_without_indels = align_len_without_indels / float(line[10])
            
            #Aligment length must be higher of 10% of the assembly transcript
            if ass_cov < 0.5:
                if ass_cov >= 0.1: trans_with_hit_list.append(line[9])
                false_pos_list.append(line[9])
                continue
            
            #Full recontructed 99%
            if ref_cov>=0.99 and ass_cov>=0.99 and indel_p<=0.01:
                if line[13] in full_rec99_dic:
                    full_rec99_dic[line[13]].append(line[9])
                else:
                    full_rec99_dic[line[13]] = [line[9]]
            #Full recontructed 95%
            if ref_cov>=0.95 and ass_cov>=0.95 and indel_p<=0.01:
                if line[13] in full_rec95_dic:
                    full_rec95_dic[line[13]].append(line[9])
                else:
                    full_rec95_dic[line[13]] = [line[9]]
            #Full recontructed 90%
            if ref_cov>=0.90 and ass_cov>=0.90 and indel_p<=0.01:
                if line[13] in full_rec90_dic:
                    full_rec90_dic[line[13]].append(line[9])
                else:
                    full_rec90_dic[line[13]] = [line[9]]
            #Full recontructed 80%
            if ref_cov>=0.80 and ass_cov>=0.80 and indel_p<=0.01:
                if line[13] in full_rec80_dic:
                    full_rec80_dic[line[13]].append(line[9])
                else:
                    full_rec80_dic[line[13]] = [line[9]]
            #Full recontructed 70%
            if ref_cov>=0.70 and ass_cov>=0.70 and indel_p<=0.01:
                if line[13] in full_rec70_dic:
                    full_rec70_dic[line[13]].append(line[9])
                else:
                    full_rec70_dic[line[13]] = [line[9]]
            
            #If the amount of indels is more of 5% is a gapped contig
            if indel_p<=0.05:
                true_pos_list.append(line[9])
                if line[9] in true_max_match_len_dic:
                    true_max_match_len_dic[line[9]] = max(true_max_match_len_dic[line[9]], int(line[0]))
                else:
                    true_max_match_len_dic[line[9]] = int(line[0])
                #If the reference coverage is 40% or more, the match is stored to quimeric and collapsed analysis
                if ref_cov>=0.40:
                    if line[9] in true_match_dic:
                        true_match_dic[line[9]].append([line[13], align_len, float(line[11]), float(line[12])])
                    else:
                        true_match_dic[line[9]] = [[line[13], align_len, float(line[11]), float(line[12])]]
            else:
                load_in_ass_gapped_dic(ass_gapped_dic, line[9], ass_cov, ass_cov_without_indels, indel_p, line[13], line[18], line[20])
                
    ass_to_ref.close()
    
    true_ass_align_len = float(sum(true_max_match_len_dic.values()))
    out_file.write("Total assembly align bps = %s\n" % true_ass_align_len)
    
    ass_len = float(commands.getoutput("infoseq -nodescription -nousa -nodatabase -noaccession -nogi -noseqversion -notype -noorganism -noheading -nopgc -noname -sequence %s | awk '{total_length=total_length+$1}END{print total_length}'" % (args.assembly_file)).split('\n')[1])
    nucleotide_especificity = true_ass_align_len / ass_len
    out_file.write("Nucleotide especificity = %s\n" % nucleotide_especificity)
    
    ref_len = float(commands.getoutput("infoseq -nodescription -nousa -nodatabase -noaccession -nogi -noseqversion -notype -noorganism -noheading -nopgc -noname -sequence %s | awk '{total_length=total_length+$1}END{print total_length}'" % (args.ref_file)).split('\n')[1])
    nucleotide_sensitivity = true_ass_align_len / ref_len
    out_file.write("Nucleotide sensitivity = %s\n" % nucleotide_sensitivity)
    
    trans_with_hit_set = set(trans_with_hit_list)
    trans_with_hit_N = float(len(trans_with_hit_set))
    out_file.write("Transcripts with hit = %s\n" % trans_with_hit_N)

    #Use set to select unique elements in different lists
    true_pos_set = set(true_pos_list)
    false_pos_set = set(false_pos_list)
    only_false_set = false_pos_set - true_pos_set
    ass_gapped_set = set(ass_gapped_dic.keys())
    only_gapped_set = ass_gapped_set - true_pos_set
    false_without_gapped_set = only_false_set - only_gapped_set

    full_rec99_set = set(full_rec99_dic.keys())
    full_rec95_set = set(full_rec95_dic.keys())
    full_rec90_set = set(full_rec90_dic.keys())
    full_rec80_set = set(full_rec80_dic.keys())
    full_rec70_set = set(full_rec70_dic.keys())

    full_rec99 = float(len(full_rec99_set)) / ref_trans_N
    full_rec95 = float(len(full_rec95_set)) / ref_trans_N
    full_rec90 = float(len(full_rec90_set)) / ref_trans_N
    full_rec80 = float(len(full_rec80_set)) / ref_trans_N
    full_rec70 = float(len(full_rec70_set)) / ref_trans_N
    
    out_file.write("Full reconstructed transcripts with 99 percent of coverage = %s\n" % full_rec99)
    out_file.write("Full reconstructed transcripts with 95 percent of coverage = %s\n" % full_rec95)
    out_file.write("Full reconstructed transcripts with 90 percent of coverage = %s\n" % full_rec90)
    out_file.write("Full reconstructed transcripts with 80 percent of coverage = %s\n" % full_rec80)
    out_file.write("Full reconstructed transcripts with 70 percent of coverage = %s\n" % full_rec70) 

    false_pos = float(len(false_without_gapped_set)) / ass_trans_N
    out_file.write("False positives = %s\n" % false_pos)
    
    gapped_contigs = float(len(only_gapped_set)) / ass_trans_N
    out_file.write("Gapped contigs = %s\n" % gapped_contigs)
    
    false_pos = open('%s_false_positives.list' % args.prefix, 'w')
    false_pos.write('%s' % '\n'.join(list(false_without_gapped_set)))
    false_pos.close()
    ass_gapped = open('%s_gapped_contigs.list' % args.prefix, 'w')
    for ass_trans in only_gapped_set:
        aux = [str(x).replace(',','-') for x in ass_gapped_dic[ass_trans]]
        ass_gapped.write('%s,%s\n' % (ass_trans, ','.join(aux)))
    ass_gapped.close()
    full_rec99 = open('%s_full_reconstructed_99.list' % args.prefix, 'w')
    for ref_trans in full_rec99_dic.keys():
        full_rec99.write('%s,%s\n' % (ref_trans, ','.join(full_rec99_dic[ref_trans])))
    full_rec99.close()
    full_rec90 = open('%s_full_reconstructed_90.list' % args.prefix, 'w')
    for ref_trans in full_rec90_dic.keys():
        full_rec90.write('%s,%s\n' % (ref_trans, ','.join(full_rec90_dic[ref_trans])))
    full_rec90.close()
    full_rec70 = open('%s_full_reconstructed_70.list' % args.prefix, 'w')
    for ref_trans in full_rec70_dic.keys():
        full_rec70.write('%s,%s\n' % (ref_trans, ','.join(full_rec70_dic[ref_trans])))
    full_rec70.close()

    full70allcontigs = []
    for full70contigs in full_rec70_dic.values():
        full70allcontigs = full70allcontigs + full70contigs
    
    #######################################################################################################################
    #QUIMERA IDENTIFICATION (MULTI GENE)
    #######################################################################################################################
    
    trans_multi_quimeric_list = []
    trans_collapsed_dic = {}
    for trans in true_match_dic:
        if len(true_match_dic[trans]) > 1:
	    max_index = map(lambda x: x[1], true_match_dic[trans]).index(max(map(lambda x: x[1], true_match_dic[trans])))
            max_match = true_match_dic[trans][max_index]
            for match in true_match_dic[trans]:
                #No view cases of alignments with the subject for the best match
                if match[0] == max_match[0]:
                    continue
                #case 1:
                #   max_match    |------------------------------------------------|
                #   match                                              |--------------------------------|  
                if match[2] > max_match[2] and match[2] < max_match[3] and match[3] > max_match[3]:
                    #If the overlap among matches is 20% or less is a quimera 
                    if max_match[3]-match[2] <= min(0.2*max_match[1], 0.2*match[1], float(100)):
                        trans_multi_quimeric_list.append(trans)
                    #If the overlap is 80% or more is a collapsed transcript  
                    if max_match[3]-match[2] >= min(0.8*max_match[1], 0.8*match[1]):
                        if trans in trans_collapsed_dic:
                            trans_collapsed_dic[trans].append(match[0])
                        else:
                            trans_collapsed_dic[trans] = [match[0]]
                            trans_collapsed_dic[trans].append(max_match[0])
                #case 2:
                #   max_match                           |------------------------------------------------|
                #   match        |--------------------------------|  
                if match[2] < max_match[2] and match[3] > max_match[2] and match[3] < max_match[3]:
                    #If the overlap among matches is 20% or less is a quimera 
                    if match[3]-max_match[2] <= min(0.2*max_match[1], 0.2*match[1], float(100)):
                        trans_multi_quimeric_list.append(trans)
                    #If the overlap is 80% or more is a collapsed transcript  
                    if match[3]-max_match[2] >= min(0.8*max_match[1], 0.8*match[1]):
                        if trans in trans_collapsed_dic:
                            trans_collapsed_dic[trans].append(match[0])
                        else:
                            trans_collapsed_dic[trans] = [match[0]]
                            trans_collapsed_dic[trans].append(max_match[0])
                #case 3:
                #   max_match    |------------------------------------------------|
                #   match                                                               |--------------------------------|  
                # or
                #   max_match                                           |------------------------------------------------|
                #   match        |--------------------------------|  
                if max_match[3] < match[2] or max_match[2] >  match[3]:
                    trans_multi_quimeric_list.append(trans)
	
    #######################################################################################################################
    #QUIMERA IDENTIFICATION (SELF GENE)                            
    #######################################################################################################################
    
    trans_self_quimeric_list = []
    #Parsing assembly to assembly psl file
    ass_to_ass = open('%s_assembly_to_assembly.psl' % args.prefix, 'r')    
    line_flag = 0
    for line in ass_to_ass:
       line_flag+=1
        #No work with first 5 lines of header file
       if line_flag>5:
           line = line.split()
           if line[9] == line[13]:
               #If align with yourself in different positions, then is a self quimera
               if line[11] != line[15] and line[12] != line[16] and float(line[0]) >= 0.2*float(line[10]):
                   trans_self_quimeric_list.append(line[9])
    
    ass_to_ass.close()
    
    #If 2 or more possible reference transcripts were reconstructed at 70% for different contigs, then, they don't are collapsed 
    trans_collapsed_dic_filter = {}
    for ass_contig in trans_collapsed_dic.keys():
        full70_list = []
        collap_flag = 0
        for ref_trans in trans_collapsed_dic[ass_contig]:
            try:
                #Not filter reference transcripts that were recontructed with 2 or more contigs at 70%
                if len(full_rec70_dic[ref_trans]) == 1:
                    full70_list=full70_list+full_rec70_dic[ref_trans]
                else:
                    collap_flag = 1
            #Not filter reference transcripts without hits on 70% reconstructed reference set
            except:
                collap_flag = 1

        #Filter contigs using conditions fixed above or if each reference transcript were recontructed with only a one different contig
        if collap_flag or (len(set(full70_list)) != len(trans_collapsed_dic[ass_contig])):
            trans_collapsed_dic_filter[ass_contig] = trans_collapsed_dic[ass_contig]

    trans_collapsed_set = set(trans_collapsed_dic_filter.keys())
    trans_multi_quimeric_set = set(trans_multi_quimeric_list)
    trans_self_quimeric_set = set(trans_self_quimeric_list)
   
    true_pos_N = float(len(true_pos_set))
    trans_collapsed = float(len(trans_collapsed_set)) / true_pos_N
    trans_multi_quimeric = float(len(trans_multi_quimeric_set)) / true_pos_N
    trans_self_quimeric = float(len(trans_self_quimeric_set)) / true_pos_N
    trans_total_quimeric = trans_multi_quimeric + trans_self_quimeric

    collap_trans = open('%s_collapsed_transcripts.list' % args.prefix, 'w')
    for ass_contig in trans_collapsed_dic_filter.keys():
        collap_trans.write('%s,%s\n' % (ass_contig, ','.join(trans_collapsed_dic_filter[ass_contig])))
    collap_trans.close()

    out_file.write("Quimeric transcripts (multiple genes) = %s\n" % trans_multi_quimeric)
    out_file.write("Quimeric transcripts (self gene) = %s\n" % trans_self_quimeric)
    out_file.write("Quimeric transcripts (total) = %s\n" % trans_total_quimeric)
    out_file.write("Collapsed contigs = %s\n" % trans_collapsed)

    #######################################################################################################################
    #FRAGMENTED REFERENCE TRANSCRIPTS
    #######################################################################################################################
    
    #Parsing reference to assembly psl file
    ref_to_ass = open('%s_reference_to_assembly.psl' % args.prefix, 'r')
    

    line_flag = 0
    ref_trans_dic = {}
    ref_ass_true_match = {}
    ref_trans_with_hit_list = []
    ref_trans_with_hit50_list = []

    for line in ref_to_ass:
        line_flag+=1
        #No work with first 5 lines of header file
        if line_flag>5:
            line = line.split()
            
            #Calculate indel and alignment lengths
            indel_N = max(float(line[5]), float(line[7]))
            align_len = float(line[0]) + float(line[1]) + float(line[2]) + indel_N
            indel_p = indel_N / align_len

            #Calculate alignment coverage over reference and assembly transcripts
            ref_cov = align_len / float(line[10])
            ass_cov = align_len / float(line[14])
            
            if line[13] in list(only_gapped_set) or line[13] in list(false_without_gapped_set): continue
            #Alignment coverage 10% or more of the reference transcript and 50% or more of the assembly contig to analize fragmented transcripts
            if ref_cov >= 0.1 and ass_cov >= 0.5 and indel_p <= 0.05:
                ref_trans_with_hit_list.append(line[9])
                if line[13] not in full70allcontigs:
                    if line[9] in ref_trans_dic:
                        ref_trans_dic[line[9]].append([line[13], align_len, float(line[11]), float(line[12]), float(line[10])])
                    else:
                        ref_trans_dic[line[9]] = [[line[13], align_len, float(line[11]), float(line[12]), float(line[10])]]
                
                #Reference coverage more of 0.5 and indels porcentage less than 1% to true match 
                if ref_cov >= 0.5 and indel_p <= 0.01:
                    ref_trans_with_hit50_list.append(line[9])
                    if line[9] in ref_ass_true_match:
                        ref_ass_true_match[line[9]].append(line[13])
                    else:
                        ref_ass_true_match[line[9]] = [line[13]]
                    
    ref_to_ass.close() 
    
    #Write list of reference transcripts names on file
    ref_with_hit50 = open('%s_reference_with_hit.list' % args.prefix, 'w')
    ref_with_hit50.write('\n'.join(list(set(ref_trans_with_hit50_list))))
    ref_with_hit50.close()
    ref_with_hit50_set = set(ref_trans_with_hit50_list)
    
    ref_trans_frag_dic = {}
    ref_trans_frag_ass_dic = {}
    for ref_trans in ref_trans_dic:
        if ref_trans in full_rec70_set: continue
        #At least 2 match for reference transcript to analize
        if len(ref_trans_dic[ref_trans]) > 1:
            matches_list = []
            #Compare all matches among them
            index = 0 
            for f_match in ref_trans_dic[ref_trans]:
                index+=1
                for match in ref_trans_dic[ref_trans][index:]:
                    #No use cases of alignments with the subject for the best match or yourself
                    if match[0] == f_match[0]:
                        continue
                    #case 1:
                    #   f_match    |------------------------------------------------|
                    #   match                                              |--------------------------------|  
                    if match[2] > f_match[2] and match[2] < f_match[3] and match[3] > f_match[3] and f_match[3]-match[2] <= min(0.2*f_match[1], 0.2*match[1], float(100)):
                        matches_list.append(match[0]) 
                        matches_list.append(f_match[0]) 
                        if ref_trans in ref_trans_frag_dic:
                            ref_trans_frag_dic[ref_trans] = [min(ref_trans_frag_dic[ref_trans][0], match[2]), max(ref_trans_frag_dic[ref_trans][1], f_match[3]), match[4]]
                        else:
                            ref_trans_frag_dic[ref_trans] = [f_match[2], match[3],  match[4]]
                        
                    #case 2:
                    #   f_match                                |------------------------------------------------|
                    #   match        |--------------------------------|  
                    if match[2] < f_match[2] and match[3] > f_match[2] and match[3] < f_match[3] and match[3]-f_match[2] <= min(0.2*f_match[1], 0.2*match[1], float(100)):
                        matches_list.append(match[0]) 
                        matches_list.append(f_match[0]) 
                        if ref_trans in ref_trans_frag_dic:
                            ref_trans_frag_dic[ref_trans] = [min(ref_trans_frag_dic[ref_trans][0], f_match[2]), max(ref_trans_frag_dic[ref_trans][1], match[3]), match[4]]
                        else:
                            ref_trans_frag_dic[ref_trans] = [match[2], f_match[3],  match[4]]
                    #case 3:
                    #   f_match    |------------------------------------------------|
                    #   match                                                               |--------------------------------|  
                    # or
                    #   f_match                                                 |------------------------------------------------|
                    #   match        |--------------------------------|  
                    if f_match[3] < match[2] or f_match[2] >  match[3]:
                        matches_list.append(match[0]) 
                        matches_list.append(f_match[0]) 
                        if ref_trans in ref_trans_frag_dic:
                            ref_trans_frag_dic[ref_trans] = [min(ref_trans_frag_dic[ref_trans][0], match[2], f_match[2]), max(ref_trans_frag_dic[ref_trans][1], match[3], f_match[3]), match[4]]
                        else:
                            ref_trans_frag_dic[ref_trans] = [min(match[2], f_match[2]), max(match[3], f_match[3]),  match[4]]

            ref_trans_frag_ass_dic[ref_trans] = list(set(matches_list))
  
    ref_trans_frag_list = []
    for ref_trans in ref_trans_frag_dic:
        ref_cov = (ref_trans_frag_dic[ref_trans][1] - ref_trans_frag_dic[ref_trans][0]) / ref_trans_frag_dic[ref_trans][2]
        #If all matches togheter covering more than 50% of the reference transcript, then it is fragmented
        if ref_trans_frag_ass_dic[ref_trans] > 1:
            if ref_cov >= 0.5:
                ref_trans_frag_list.append(ref_trans)
        else:
                ref_trans_frag_list.append(ref_trans)


    ref_trans_with_hit_set = set(ref_trans_with_hit_list)
    ref_trans_frag = float(len(ref_trans_frag_list)) / float(len(ref_trans_with_hit_set))
    out_file.write("Fragmented reference transcripts = %s\n" % ref_trans_frag)

    frag_trans = open('%s_fragmented_transcripts.list' % args.prefix, 'w')
    for ref_trans in ref_trans_frag_list:
        frag_trans.write('%s,%s\n' % (ref_trans, ','.join(ref_trans_frag_ass_dic[ref_trans])))
    frag_trans.close()

    #######################################################################################################################
    #ISOFORMS IDENTIFICATION
    #######################################################################################################################

    #Only oases and trinity provide isoforms identification
    if args.a_soft in ['oases', 'trinity', 'bridger']:
        
        #Parsing headers of assembled transcripts saving gene name and transcript amount for each gene
        #all_ass_trans = commands.getoutput("fgrep '>' %s | cut -d ' ' -f 1 | sed 's/>//g'" % (args.assembly_file)).split('\n')
        all_ass_gene_dic = {}
        for ass_trans in all_ass_trans:
            if args.a_soft == "trinity":
                ass_gene = "_".join(ass_trans.split('_')[0:4])
            if args.a_soft == "oases":
                ass_gene = "_".join(ass_trans.split('_')[0:2])
            if args.a_soft == "bridger":
                ass_gene = ass_trans.split('_')[0]
            if ass_gene in all_ass_gene_dic:
                all_ass_gene_dic[ass_gene]+=1
            else:
                all_ass_gene_dic[ass_gene]=1

        #Parsing gene transcripts reference information file
        iso_file = open(args.gene_trans_file, 'r')
        gene_trans_dic = {}
        for line in iso_file:
            gene,trans = line.split(',')
            if gene in gene_trans_dic:
                gene_trans_dic[gene].append(trans[:-1])
            else:
                gene_trans_dic[gene] = [trans[:-1]]
        iso_file.close()        
        
        one_gene_correct_iso_list = []
        one_gene_more_iso_list = []
        one_gene_less_iso_list = []
        more_one_gene_correct_iso_list = []
        more_one_gene_more_iso_list = []
        more_one_gene_less_iso_list = []
        ref_gene_ass = []
 
        for gene in gene_trans_dic:
            ass_trans_list = []
            for trans in gene_trans_dic[gene]:
                #Saving the name of assembled transcripts that align with each transcript of this gene
                if trans in ref_ass_true_match:
                    ass_trans_list = ass_trans_list + ref_ass_true_match[trans]
            #If at least one assembled transcript align
            if ass_trans_list:
                #Saving the name of the reference gene that a assembled gene
                ref_gene_ass.append(gene)
                ass_trans_set = set(ass_trans_list)
                
                #Parsing names of assembled transcripts to identify genes and isoforms
                ass_gene_dic = {}
                for ass_trans in ass_trans_set:
                    if args.a_soft == "trinity":
                        ass_gene = "_".join(ass_trans.split('_')[0:4])
                    if args.a_soft == "oases":
                        ass_gene = "_".join(ass_trans.split('_')[0:2])
                    if args.a_soft == "bridger":
                        ass_gene = ass_trans.split('_')[0]
                    #Saving the amount of predicted isoforms for each assembled 'gene' 
                    if ass_gene in ass_gene_dic:
                        ass_gene_dic[ass_gene]+=1
                    else:
                        ass_gene_dic[ass_gene]=1
                #If the transcripts for this gene align only with assembled transcripts of one predicted 'gene'
                if len(ass_gene_dic) == 1:
                    #If the amount of predicted isoforms is the same of the reference then it is the correct situation
                    if len(gene_trans_dic[gene]) == ass_gene_dic.values()[0] and len(gene_trans_dic[gene]) == all_ass_gene_dic[ass_gene_dic.keys()[0]]:
                        one_gene_correct_iso_list.append(gene)
                    #Predicted more isoforms that the reference 
                    elif len(gene_trans_dic[gene]) < ass_gene_dic.values()[0] or len(gene_trans_dic[gene]) < all_ass_gene_dic[ass_gene_dic.keys()[0]]:
                        one_gene_more_iso_list.append(gene)
                    #Predicted less isoforms that the reference 
                    elif len(gene_trans_dic[gene]) > ass_gene_dic.values()[0] and len(gene_trans_dic[gene]) > all_ass_gene_dic[ass_gene_dic.keys()[0]]:
                        one_gene_less_iso_list.append(gene)
                #If the transcripts for this gene align with assembled transcripts for several predicted 'genes'
                else:
                    if len(gene_trans_dic[gene]) == sum(ass_gene_dic.values()) and len(gene_trans_dic[gene]) == sum([all_ass_gene_dic[gene_ass] for gene_ass in ass_gene_dic.keys()]):
                        more_one_gene_correct_iso_list.append(gene)
                    elif len(gene_trans_dic[gene]) < sum(ass_gene_dic.values()) or len(gene_trans_dic[gene]) < sum([all_ass_gene_dic[gene_ass] for gene_ass in ass_gene_dic.keys()]):
                        more_one_gene_more_iso_list.append(gene)
                    elif len(gene_trans_dic[gene]) > sum(ass_gene_dic.values()) and len(gene_trans_dic[gene]) > sum([all_ass_gene_dic[gene_ass] for gene_ass in ass_gene_dic.keys()]):
                        more_one_gene_less_iso_list.append(gene)
            
        ref_gene_ass_N = float(len(set(ref_gene_ass)))
        ref_gene_N = len(gene_trans_dic)
        ass_refgene_relation = ref_gene_ass_N / ref_gene_N
        one_gene_correct_iso = float(len(one_gene_correct_iso_list)) / ref_gene_ass_N
        one_gene_more_iso = float(len(one_gene_more_iso_list)) / ref_gene_ass_N
        one_gene_less_iso = float(len(one_gene_less_iso_list)) / ref_gene_ass_N
        more_one_gene_correct_iso = float(len(more_one_gene_correct_iso_list)) / ref_gene_ass_N
        more_one_gene_more_iso = float(len(more_one_gene_more_iso_list)) / ref_gene_ass_N     
        more_one_gene_less_iso = float(len(more_one_gene_less_iso_list)) / ref_gene_ass_N     

        out_file.write("Number of reference genes = %s\n" % ref_gene_N)
        out_file.write("Number of assembled genes = %s\n" % ref_gene_ass_N)
        out_file.write("Assembled genes / Reference genes = %s\n" % ass_refgene_relation)
        out_file.write("Reference genes with correct predicted amount of transcripts (1 gene) = %s\n" % one_gene_correct_iso)
        out_file.write("Reference genes with more predicted amount of transcripts (1 gene) = %s\n" % one_gene_more_iso)
        out_file.write("Reference genes with less predicted amount of transcripts (1 gene) = %s\n" % one_gene_less_iso)
        out_file.write("Reference genes with correct predicted amount of transcripts (2 or more genes) = %s\n" % more_one_gene_correct_iso)
        out_file.write("Reference genes with more predicted amount of transcripts (2 or more genes) = %s\n" % more_one_gene_more_iso)
        out_file.write("Reference genes with less predicted amount of transcripts (2 or more genes) = %s\n" % more_one_gene_less_iso)

    #######################################################################################################################
    #TRANSCRIPT ASSEMBLED WITH LOW EXPRESSION
    #######################################################################################################################
    
    exp_file = open(args.exp_file, 'r')
    #exp_file = reference_transcript,expression_value

    exp_dict = {}
    for line in exp_file:
        line = line[:-1].split(',')
        exp_dict[line[0]] = float(line[1])
    exp_file.close()
    #If expression value is low of 10 percentile of expression distribution values
    percentile10_value = float(np.percentile(exp_dict.values(), 10))
    ref_trans_p10_exp = [trans for trans in exp_dict.keys() if exp_dict[trans] <= percentile10_value]
    
    cutoff_value = float(10)
    ref_trans_ctv_exp = [trans for trans in exp_dict.keys() if exp_dict[trans] <= cutoff_value]

    #If the reference transcript with low expression was assembled
    ref_trans_p10_exp_ass = [trans for trans in ref_trans_p10_exp if trans in full_rec90_set]
    ref_trans_ctv_exp_ass = [trans for trans in ref_trans_ctv_exp if trans in full_rec90_set]
    
    ref_trans_p10_exp_ass_N = float(len(ref_trans_p10_exp_ass))
    ref_trans_p10_exp_ass = ref_trans_p10_exp_ass_N / float(len(ref_trans_p10_exp))
    out_file.write("Reference transcripts with expression in the 10 percentile assembled = %s\n" % ref_trans_p10_exp_ass)
    
    ref_trans_ctv_exp_ass_N = float(len(ref_trans_ctv_exp_ass))
    try:
        ref_trans_ctv_exp_ass = ref_trans_ctv_exp_ass_N / float(len(ref_trans_ctv_exp))
    except:
        ref_trans_ctv_exp_ass = float(0)
    out_file.write("Reference transcripts with expression less than %s assembled = %s\n" % (cutoff_value, ref_trans_ctv_exp_ass))

    out_file.close()

    #######################################################################################################################
    #REMOVE .psl FILES
    #######################################################################################################################
    #os.system("rm *.psl")
