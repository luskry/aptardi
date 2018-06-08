#!/bin/bash python3

import re
import os
import datetime
from bisect import bisect
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import fbeta_score, make_scorer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import NearMiss
import pickle

temp_dir = '/Users/ryan/Desktop/ml_prototype/'
gold_standard_rna_seq_sorted_bam_file = temp_dir+'sorted_uhr_exp2_auto_ensembl_v91_genome_index_default_hisat2.bam'
gold_standard_polya_seq_sorted_bam_file = temp_dir+'sorted_all_uhr_polya_trimmed_ensembl_v91_genome_index_default_hisat2.bam'
gold_standard_gtf_file_from_rna_seq_data = temp_dir+'uhr_exp2_auto_ensembl_v91_genome_index_default_GRCh38.v91.txt'
chrom_sizes_file = temp_dir+'hg38_self_built.chrom.sizes'
fasta_file = temp_dir+'genome_edited_headers.fa'
sequence_bed_file = 'modified_bed_with_sequence.bed'
UTR_bed_file = 'UTR_seq.bed' 
All_Wig_files = ['gold_standard_rna_seq_coverage.bedgraph']
All_Wig_files_polya = ['gold_standard_polya_seq_coverage.bedgraph']

def time_now():
    
    curr_time = datetime.datetime.now()
    return curr_time.strftime('%c')   


def gtf_to_bed(gtf_file, tmp_dir):
    
    
    def bed_line(estart, eend, field, nline, exon_coverage, trans_coverage):
        
        allids = {}
        estp = estart[0] - 1
        eedp = eend[-1]
        geneid = re.findall(r'gene_name \"([\w\.]+)\"', field[8])
        if len(geneid) == 0:
            geneid = re.findall(r'ref_gene_name \"([\w\.]+)\"', field[8])
        if len(geneid) == 0:
            geneid = re.findall(r'gene_id \"([\w\.]+)\"', field[8])   
        if len(geneid) == 0:
            geneid = 'NA'
        else:
            geneid = geneid[0]
        transid = re.findall(r'reference_id \"([\w\.]+)\"', field[8])
        if len(transid) == 0:
            transid = re.findall(r'transcript_id \"([\w\.]+)\"', field[8])
        if len(transid) == 0:
            transid = 'Trans_' + str(nline)
        else:
            transid = transid[0]
        if transid in allids.keys():
            transid2 = transid + '_DUP' + str(allids[transid])
            allids[transid] = allids[transid] + 1
            transid = transid2
        else:
            allids[transid] = 1
        if not field[0].startswith('chr'):
            field[0] = 'chr' + field[0]
        line = [field[0], str(estp), str(eedp), '|'.join([transid, geneid, field[0], str(estp), str(eedp), field[6]]), ','.join(exon_coverage), field[6], str(estp), str(eedp), str(trans_coverage), str(len(estart))]
        seglen = [eend[i] - estart[i] + 1 for i in range(len(estart))]
        segstart = [estart[i] - estart[0] for i in range(len(estart))]
        strl = str(seglen[0])
        for i in range(1, len(seglen)):
            strl += ',' + str(seglen[i])
        strs = str(segstart[0])
        for i in range(1, len(segstart)):
            strs +=',' + str(segstart[i])
        line.extend([strl, strs])
        
        return line
    
    
    extracted_bed_line = []   
    estart = []
    eend = []
    exon_coverage = []
    trans_coverage = []
    nline = 0
    prevfield = []
    prevtransid = ''
    with open(gtf_file, 'r') as f:
        for lines in f:
            if lines.strip().split('\t')[0] == '\n' or lines[0] == '#':
                continue
            else:
                field = lines.strip().split('\t')
                field[8:] = [(' ').join(field[8:])]
                nline = nline + 1
                if field[1] not in ['Cufflinks', 'StringTie']:
                    print('Warning: The GTF file may not have been produced by Cufflinks or StringTie and therefore may give erroneous results')
                    pass
                if field[2] != 'exon' and field[2] != 'transcript':
                    continue
                transid = re.findall(r'transcript_id \"([\w\.]+)\"', field[8])
                if len(transid) > 0:
                    transid = transid[0]
                else:
                    transid = ''
                if field[2] == 'transcript' or (prevtransid != '' and transid != '' and transid != prevtransid):
                    if len(estart) != 0:
                        extracted_bed_line.append(bed_line(estart, eend, prevfield, nline, exon_coverage, trans_coverage))
                    estart = []
                    eend = []
                    exon_coverage = []
                    trans_coverage = []
                prevfield = field
                prevtransid = transid
                if field[2] == 'exon':
                    est = int(field[3])
                    eed = int(field[4])
                    exon_cov = re.findall(r'cov \"([\d\.]+)\"', field[8])
                    if len(exon_cov) == 0:
                        exon_cov = ['NA']
                    estart += [est]
                    eend += [eed]
                    exon_coverage.extend(exon_cov)
                if field[2] == 'transcript':
                    trans_cov = re.findall(r'cov \"([\d\.]+)\"', field[8])
                if len(trans_cov) == 0:
                    trans_coverage = '255,0,0'
                else:
                    trans_coverage = trans_cov[0]
        if len(estart) != 0:
            extracted_bed_line.append(bed_line(estart, eend, prevfield, nline, exon_coverage, trans_coverage))
                
    return extracted_bed_line
    
    
bed_list = gtf_to_bed(gold_standard_gtf_file_from_rna_seq_data, temp_dir)  

    
def bam_to_bedgraph(rna_seq_bam_path, polya_seq_bam_path, tmp_dir):
      
    os.system('bedtools genomecov -bg -strand + -ibam %s > %s' % (rna_seq_bam_path, tmp_dir+'rna_seq_coverage_plus_strand.bedgraph'))
    os.system('bedtools genomecov -bg -strand - -ibam %s > %s' % (rna_seq_bam_path, tmp_dir+'rna_seq_coverage_minus_strand.bedgraph'))
    os.system('bedtools unionbedg -i %s %s > %s' % (tmp_dir+'rna_seq_coverage_plus_strand.bedgraph', tmp_dir+'rna_seq_coverage_minus_strand.bedgraph', tmp_dir+'gold_standard_rna_seq_coverage.bedgraph'))
        
    try:
        os.remove(tmp_dir+'rna_seq_coverage_plus_strand.bedgraph')
    except OSError:
        pass
    
    try:
        os.remove(tmp_dir+'rna_seq_coverage_minus_strand.bedgraph')
    except OSError:
        pass
      
    os.system('bedtools genomecov -bg -strand + -5 -ibam %s > %s' % (polya_seq_bam_path, tmp_dir+'polya_seq_coverage_plus_strand.bedgraph'))
    os.system('bedtools genomecov -bg -strand - -5 -ibam %s > %s' % (polya_seq_bam_path, tmp_dir+'polya_seq_coverage_minus_strand.bedgraph'))
    os.system('bedtools unionbedg -i %s %s > %s' % (tmp_dir+'polya_seq_coverage_plus_strand.bedgraph', tmp_dir+'polya_seq_coverage_minus_strand.bedgraph', tmp_dir+'gold_standard_polya_seq_coverage.bedgraph'))
    
    try:
        os.remove(tmp_dir+'polya_seq_coverage_plus_strand.bedgraph')
    except OSError:
        pass
    
    try:
        os.remove(tmp_dir+'polya_seq_coverage_minus_strand.bedgraph')
    except OSError:
        pass
    
    
bam_to_bedgraph(gold_standard_rna_seq_sorted_bam_file, gold_standard_polya_seq_sorted_bam_file, temp_dir)


def get_chrom_sizes(chrom_file):
    
    chrom_sizes_dict = {}
    
    if chrom_file.endswith('.gz'):
        os.system('gunzip -f %s' % (chrom_file))
        chrom_file = chrom_file[:-3]
        
    with open(chrom_file, 'r') as f:
        for line in f:
            fields = line.strip('\n').split('\t')
            if line.strip().split('\t')[0] == '\n' or line[0] == '#':
                pass
            elif '_' not in fields[0]:
                chrom = fields[0]
                chrom_size = fields[1]
                if not chrom.startswith('chr'):
                    chrom = 'chr' + chrom
            chrom_sizes_dict[chrom] = [chrom_size]
            
    return chrom_sizes_dict


chrom_dict = get_chrom_sizes(chrom_sizes_file)


def annot_from_bed(bed_list_raw, chrom_sizes_dict, tmp_dir):
    
    utr_size = 2000
    output_write = open(tmp_dir+'modified_annot.bed', 'w')
    
    for fields in bed_list_raw:
        if '_' not in fields[0]:
            strand = fields[5]
            tx_start = fields[1]
            tx_end = fields[2]  
            abs_exon_starts = [int(tx_start)+int(i) for i in fields[11].strip(',').split(',')]
            exon_sizes = [int(i) for i in fields[10].strip(',').split(',')]
        
            if strand == '+' and int(fields[10].strip(',').split(',')[-1]) < utr_size:
                extended_UTR = str(int(tx_end) + (utr_size - int(fields[10].strip(',').split(',')[-1]))) 
                if fields[0] in chrom_sizes_dict:
                    total_chrom_length = chrom_sizes_dict[fields[0]][0]
                    if int(total_chrom_length) < int(extended_UTR):
                        extended_UTR = total_chrom_length
                tx_end = extended_UTR
                exon_sizes[-1] = int(tx_end)-abs_exon_starts[-1]
                        
            elif strand == '-' and int(fields[10].strip(',').split(',')[0]) < utr_size:
                extended_UTR = str(int(tx_start) - (utr_size - int(fields[10].strip(',').split(',')[0])))
                if int(extended_UTR) < 0:
                    extended_UTR = '0'
                tx_start = extended_UTR
                abs_exon_starts[0] = int(extended_UTR)
                exon_sizes[0] = int(fields[1])-int(tx_start)+exon_sizes[0]
                
            else:
                continue
            
            rel_exon_starts = (',').join([str(int(i)-int(tx_start)) for i in abs_exon_starts])
            exon_sizes = (',').join(str(i) for i in exon_sizes)
            write_line = [fields[0], tx_start, tx_end, fields[3], '0', strand, fields[6], fields[7], '255, 0, 0', fields[9], exon_sizes, rel_exon_starts]
            output_write.writelines('\t'.join(write_line) + '\n')
     
    output_write.close()


annot_from_bed(bed_list, chrom_dict, temp_dir)


def sequence_from_fasta(local_fasta_file, tmp_dir):
    
    if local_fasta_file.endswith('.gz'):
        os.system('gunzip -f %s' % (local_fasta_file))
        local_fasta_file = local_fasta_file[:-3]
    
    os.system('bedtools getfasta -fi %s -bed %s -name -bedOut -s > %s' % (local_fasta_file, tmp_dir+'modified_annot.bed', tmp_dir+'modified_bed_with_sequence.bed'))

    try:
        os.remove(tmp_dir+'modified_annot.bed')
    except OSError:
        pass


sequence_from_fasta(fasta_file, temp_dir)


def extract_mod_UTR_seq(sequence_bed, tmp_dir):
    
    utr_size = 2000
    output_write = open(tmp_dir+'UTR_seq.bed', 'w')
    
    with open(tmp_dir+sequence_bed, 'r') as f:
        for line in f:
            fields = line.strip('\n').split('\t')
            strand = fields[5]
            if strand == '+' and int(fields[10].strip(',').split(',')[-1]) >= utr_size:
                all_seq = fields[-1]
                seq = all_seq[-utr_size:]
                tx_end = str(int(fields[1]) + utr_size)
                fields[10] = fields[10].strip(',').split(',')[-1]
                fields[11] = fields[11].strip(',').split(',')[-1]
            elif strand == '+' and int(fields[10].strip(',').split(',')[-1]) < utr_size:
                seq = fields[-1]    
                tx_end = str(int(fields[1]) + len(seq))
                fields[10] = fields[10].strip(',').split(',')[-1]
                fields[11] = fields[11].strip(',').split(',')[-1]
            elif strand == '-' and int(fields[10].strip(',').split(',')[0]) >= utr_size:
                all_seq = fields[-1]
                seq = all_seq[-utr_size:]
                tx_end = str(int(fields[1]) + utr_size)
                fields[10] = fields[10].strip(',').split(',')[0]
                fields[11] = fields[11].strip(',').split(',')[0]
            elif strand == '-' and int(fields[10].strip(',').split(',')[0]) < utr_size:
                seq = fields[-1] 
                tx_end = str(int(fields[1]) + len(seq))
                fields[10] = fields[10].strip(',').split(',')[0]
                fields[11] = fields[11].strip(',').split(',')[0]
                
            fields[-1] = seq
            fields[2] = tx_end
            fields[9] = '1'
            output_write.writelines('\t'.join(fields) + '\n')
         
    output_write.close()
    
    try:
        os.remove(tmp_dir+sequence_bed)
    except OSError:
        pass

    
extract_mod_UTR_seq(sequence_bed_file, temp_dir)


def subtract_modified_UTR_overlap(UTR_bed_path, tmp_dir):
    
    os.system('sort -k1,1 -k2,2n -k3,3n -k6,6 -u %s > %s' % (tmp_dir+UTR_bed_path, tmp_dir+'UTR_duplicate_subtracted.bed'))
    
    try:
        os.remove(tmp_dir+UTR_bed_path)
    except OSError:
        pass


subtract_modified_UTR_overlap(UTR_bed_file, temp_dir)


def dinuc_seq_freq(UTR_bed_path, tmp_dir):
    
    UTR_dict = {}
    base_combos = ['AA', 'AC','AG','AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG','GT', 'TA', 'TC', 'TG', 'TT']
    upstream_region_one = 70
    upstream_region_two = 70
    downstream_region_one = 30
    downstream_region_two = 30
    
    with open(tmp_dir+UTR_bed_path, 'r') as f:
        for line in f:
            fields = line.strip('\n').split('\t')
            UTR_dict[fields[3]] = fields[0:3]+[fields[5]]
            bases = fields[-1]
            for i in range(len(bases)):
                if i >= upstream_region_one+upstream_region_two and i < (len(bases) - (downstream_region_one+downstream_region_two)): 
                    upstream_region_one_probe = range(i - (upstream_region_one+upstream_region_two), i - (upstream_region_two))
                    upstream_region_two_probe = range(i - upstream_region_two, i)
                    downstream_region_one_probe = range(i, i + downstream_region_one)
                    downstream_region_two_probe = range(i + downstream_region_one, i + (downstream_region_one+downstream_region_two))
                    seq_list_upstream_region_one = []
                    seq_list_upstream_region_two = []
                    seq_list_downstream_region_one = []
                    seq_list_downstream_region_two = []
                    for nuc in upstream_region_one_probe:
                        first_base = bases[nuc]
                        second_base = bases[nuc + 1]
                        dinuc = first_base+second_base
                        seq_list_upstream_region_one.append(dinuc)
                        seq_list_upstream_region_one_final = [[x, seq_list_upstream_region_one.count(x)] for x in set(seq_list_upstream_region_one)]
                        for combo in base_combos:
                            if combo not in [x[0] for x in seq_list_upstream_region_one_final]:
                                seq_list_upstream_region_one_final.append([combo, 0])
                        seq_list_upstream_region_one_final = sorted(seq_list_upstream_region_one_final, key = lambda x: x[0])
                        dinuc_sum = sum(x[1] for x in seq_list_upstream_region_one_final)
                        dinuc_freq = [x[1]/dinuc_sum for x in seq_list_upstream_region_one_final]
                        seq_list_upstream_region_one_final = [x for x in seq_list_upstream_region_one_final if 'N' not in x[0]]
                        for sub in range(len(seq_list_upstream_region_one_final)):
                            seq_list_upstream_region_one_final[sub][1] = dinuc_freq[sub]
                    for nuc in upstream_region_two_probe:
                        first_base = bases[nuc]
                        second_base = bases[nuc + 1]
                        dinuc = first_base+second_base
                        seq_list_upstream_region_two.append(dinuc)
                        seq_list_upstream_region_two_final = [[x, seq_list_upstream_region_two.count(x)] for x in set(seq_list_upstream_region_two)]
                        for combo in base_combos:
                            if combo not in [x[0] for x in seq_list_upstream_region_two_final]:
                                seq_list_upstream_region_two_final.append([combo, 0])
                        seq_list_upstream_region_two_final = sorted(seq_list_upstream_region_two_final, key = lambda x: x[0])
                        dinuc_sum = sum(x[1] for x in seq_list_upstream_region_two_final)
                        dinuc_freq = [x[1]/dinuc_sum for x in seq_list_upstream_region_two_final]
                        seq_list_upstream_region_two_final = [x for x in seq_list_upstream_region_two_final if 'N' not in x[0]]
                        for sub in range(len(seq_list_upstream_region_two_final)):
                            seq_list_upstream_region_two_final[sub][1] = dinuc_freq[sub]
                    for nuc in downstream_region_one_probe:
                        first_base = bases[nuc]
                        second_base = bases[nuc + 1]
                        dinuc = first_base+second_base
                        seq_list_downstream_region_one.append(dinuc)
                        seq_list_downstream_region_one_final = [[x, seq_list_downstream_region_one.count(x)] for x in set(seq_list_downstream_region_one)]
                        for combo in base_combos:
                            if combo not in [x[0] for x in seq_list_downstream_region_one_final]:
                                seq_list_downstream_region_one_final.append([combo, 0])
                        seq_list_downstream_region_one_final = sorted(seq_list_downstream_region_one_final, key = lambda x: x[0])
                        dinuc_sum = sum(x[1] for x in seq_list_downstream_region_one_final)
                        dinuc_freq = [x[1]/dinuc_sum for x in seq_list_downstream_region_one_final]
                        seq_list_downstream_region_one_final = [x for x in seq_list_downstream_region_one_final if 'N' not in x[0]]
                        for sub in range(len(seq_list_downstream_region_one_final)):
                            seq_list_downstream_region_one_final[sub][1] = dinuc_freq[sub]
                    for nuc in downstream_region_two_probe:
                        first_base = bases[nuc]
                        second_base = bases[nuc + 1]
                        dinuc = first_base+second_base
                        seq_list_downstream_region_two.append(dinuc)
                        seq_list_downstream_region_two_final = [[x, seq_list_downstream_region_two.count(x)] for x in set(seq_list_downstream_region_two)]
                        for combo in base_combos:
                            if combo not in [x[0] for x in seq_list_downstream_region_two_final]:
                                seq_list_downstream_region_two_final.append([combo, 0])
                        seq_list_downstream_region_two_final = sorted(seq_list_downstream_region_two_final, key = lambda x: x[0])
                        dinuc_sum = sum(x[1] for x in seq_list_downstream_region_two_final)
                        dinuc_freq = [x[1]/dinuc_sum for x in seq_list_downstream_region_two_final]
                        seq_list_downstream_region_two_final = [x for x in seq_list_downstream_region_two_final if 'N' not in x[0]]
                        for sub in range(len(seq_list_downstream_region_two_final)):
                            seq_list_downstream_region_two_final[sub][1] = dinuc_freq[sub]
                    UTR_dict[fields[3]].append(seq_list_upstream_region_one_final+seq_list_upstream_region_two_final+seq_list_downstream_region_one_final+seq_list_downstream_region_two_final)
    
    UTR_dict_final = {}
    for event in UTR_dict:
        new_list = []
        for base in range(4, len( UTR_dict[event])):
            new_list.append(np.array([j for i in UTR_dict[event][base] for j in i][1::2]))
            UTR_dict_final[event] = new_list
        
    return UTR_dict_final


seq_feat = dinuc_seq_freq('UTR_duplicate_subtracted.bed', temp_dir)


def strong_signal_variant_indicator(UTR_bed_path, tmp_dir):
    
    UTR_dict = {}
    signal_variants = ['AATAAA', 'ATTAAA']
    upstream_region_one = 70
    upstream_region_two = 70
    downstream_region_one = 30
    downstream_region_two = 30
    
    with open(tmp_dir+UTR_bed_path, 'r') as f:
        for line in f:
            fields = line.strip('\n').split('\t')
            bases = fields[-1]
            i_line = -1
            signal_variant_indicator = np.zeros(len(bases) - (upstream_region_one + upstream_region_two + downstream_region_one + downstream_region_two))
            for i in range(len(bases)):
                if i >= (upstream_region_one+upstream_region_two) and i < (len(bases) - (downstream_region_one+downstream_region_two)): 
                    probe = range(i - 30, i - 5) #window = 36th-7th base upstream base of interest not including base of interest
                    hex_seq_list = []
                    i_line += 1
                    for nuc in probe:
                        hexamer = bases[nuc-6:nuc]
                        hex_seq_list.append(hexamer)
                        if any(i in signal_variants for i in hex_seq_list):
                            signal_variant_indicator[i_line] = 1
                        else:
                            signal_variant_indicator[i_line] = 0
            UTR_dict[fields[3]] = signal_variant_indicator
            
    return UTR_dict


ssv = strong_signal_variant_indicator('UTR_duplicate_subtracted.bed', temp_dir)


def weak_signal_variant_indicator(UTR_bed_path, tmp_dir):
    
    UTR_dict = {}
    signal_variants = ['AAGAAA', 'AAAAAG', 'AATACA', 'TATAAA', 'ACTAAA', 'AGTAAA', 'GATAAA', 'AATATA', 'CATAAA', 'AATAGA']
    upstream_region_one = 70
    upstream_region_two = 70
    downstream_region_one = 30
    downstream_region_two = 30
    
    with open(tmp_dir+UTR_bed_path, 'r') as f:
        for line in f:
            fields = line.strip('\n').split('\t')
            bases = fields[-1]
            i_line = -1
            signal_variant_indicator = np.zeros(len(bases) - (upstream_region_one + upstream_region_two + downstream_region_one + downstream_region_two))
            for i in range(len(bases)):
                if i >= (upstream_region_one+upstream_region_two) and i < (len(bases) - (downstream_region_one+downstream_region_two)): 
                    probe = range(i - 30, i - 5) #window = 36th-7th base upstream base of interest not including base of interest
                    hex_seq_list = []
                    i_line += 1
                    for nuc in probe:
                        hexamer = bases[nuc-6:nuc]
                        hex_seq_list.append(hexamer)
                        if any(i in signal_variants for i in hex_seq_list):
                            signal_variant_indicator[i_line] = 1
                        else:
                            signal_variant_indicator[i_line] = 0
            UTR_dict[fields[3]] = signal_variant_indicator
            
    return UTR_dict


wsv = weak_signal_variant_indicator('UTR_duplicate_subtracted.bed', temp_dir)


def ca_element_indicator(UTR_bed_path, tmp_dir):
    
    UTR_dict = {}
    upstream_region_one = 70
    upstream_region_two = 70
    downstream_region_one = 30
    downstream_region_two = 30
    
    with open(tmp_dir+UTR_bed_path, 'r') as f:
        for line in f:
            fields = line.strip('\n').split('\t')
            bases = fields[-1]
            i_line = -1
            ca_element_indicator = np.zeros(len(bases) - (upstream_region_one + upstream_region_two + downstream_region_one + downstream_region_two))
            for i in range(len(bases)):
                if i >= (upstream_region_one+upstream_region_two) and i < (len(bases) - (downstream_region_one+downstream_region_two)): 
                    i_line += 1
                    first_base = bases[i - 2]
                    second_base = bases[i - 1]
                    dinuc = first_base+second_base
                    if dinuc == 'CA':
                        ca_element_indicator[i_line] = 1
                    else:
                        ca_element_indicator[i_line] = 0
            UTR_dict[fields[3]] = ca_element_indicator
            
    return UTR_dict


ca = ca_element_indicator('UTR_duplicate_subtracted.bed', temp_dir)


def load_polya_seq_coverage(All_Wig_files_polya, UTR_Annotation_file, tmp_dir):
    
    min_cov = 10
    UTR_events_dict = {}
    All_Samples_Total_depth = []
    All_samples_extracted_3UTR_coverage_dict = {}
    All_samples_extracted_3UTR_coverage_dict_minus = {}
    upstream_region_one = 70
    upstream_region_two = 70
    downstream_region_one = 30
    downstream_region_two = 30
    
    with open(tmp_dir+UTR_Annotation_file,'r') as f: 
        for line in f:
            fields = line.strip('\n').split('\t')
            curr_chr = fields[0]
            if not curr_chr.startswith('chr'):
                curr_chr = 'chr' + curr_chr
            region_start = fields[1]
            region_end = fields[2]
            UTR_pos = '%s:%s-%s' % (curr_chr, region_start, region_end)
            UTR_events_dict[fields[3]] = [fields[0], region_start, region_end, fields[5], UTR_pos]
       
    for curr_wig_file in All_Wig_files_polya:
        
        cur_sample_total_depth = 0
        num_line = 0
        curr_sample_All_chroms_coverage_dict = {}
        curr_sample_All_chroms_coverage_dict_minus = {}
        
        with open(tmp_dir+curr_wig_file, 'r') as f:
            for line in f:
                if '#' not in line:
                    fields = line.strip('\n').split('\t')
                    chrom_name = fields[0]
                    if not chrom_name.startswith('chr'):
                          chrom_name = 'chr' + chrom_name
                    region_start = int(fields[1])
                    region_end = int(fields[2])
                    cur_sample_total_depth += int(float(fields[-2])) * (region_end - region_start) + int(float(fields[-1])) * (region_end - region_start)
                    
                    if chrom_name not in curr_sample_All_chroms_coverage_dict:
                        curr_sample_All_chroms_coverage_dict[chrom_name] = [[0],[0]]
                    if region_start > curr_sample_All_chroms_coverage_dict[chrom_name][0][-1]:
                        curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_start)
                        curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)
                    curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_end)
                    curr_sample_All_chroms_coverage_dict[chrom_name][1].append(int(fields[-2]))
            
                    if chrom_name not in curr_sample_All_chroms_coverage_dict_minus:
                        curr_sample_All_chroms_coverage_dict_minus[chrom_name] = [[0],[0]]
                    if region_start > curr_sample_All_chroms_coverage_dict_minus[chrom_name][0][-1]:
                        curr_sample_All_chroms_coverage_dict_minus[chrom_name][0].append(region_start)
                        curr_sample_All_chroms_coverage_dict_minus[chrom_name][1].append(0)
                    curr_sample_All_chroms_coverage_dict_minus[chrom_name][0].append(region_end)
                    curr_sample_All_chroms_coverage_dict_minus[chrom_name][1].append(int(float(fields[-1])))
                num_line += 1
        All_Samples_Total_depth.append(cur_sample_total_depth)
        curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)
        curr_sample_All_chroms_coverage_dict_minus[chrom_name][1].append(0)
        
        for curr_3UTR_event_id in UTR_events_dict:
            curr_3UTR_structure = UTR_events_dict[curr_3UTR_event_id]
            curr_chr = curr_3UTR_structure[0]
            
            if UTR_events_dict[curr_3UTR_event_id][-2] == '+':
                if curr_chr in curr_sample_All_chroms_coverage_dict:
                    curr_chr_coverage = curr_sample_All_chroms_coverage_dict[curr_chr]
                    region_start = int(curr_3UTR_structure[1])
                    region_end = int(curr_3UTR_structure[2])
                    left_region_index = bisect(curr_chr_coverage[0],region_start)
                    right_region_index = bisect(curr_chr_coverage[0],region_end)
                    extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]
                    extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                    extracted_3UTR_region.insert(0,region_start)
                    extracted_3UTR_region.append(region_end)
                    
                    if curr_3UTR_event_id not in All_samples_extracted_3UTR_coverage_dict:
                        All_samples_extracted_3UTR_coverage_dict[curr_3UTR_event_id] = []
                    All_samples_extracted_3UTR_coverage_dict[curr_3UTR_event_id].append([extracted_coverage,extracted_3UTR_region])
           
            elif UTR_events_dict[curr_3UTR_event_id][-2] == '-':
                if curr_chr in curr_sample_All_chroms_coverage_dict_minus:
                    curr_chr_coverage = curr_sample_All_chroms_coverage_dict_minus[curr_chr]
                    region_start = int(curr_3UTR_structure[1])
                    region_end = int(curr_3UTR_structure[2])
                    left_region_index = bisect(curr_chr_coverage[0],region_start)
                    right_region_index = bisect(curr_chr_coverage[0],region_end)
                    extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]
                    extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                    extracted_3UTR_region.insert(0,region_start)
                    extracted_3UTR_region.append(region_end)
                    
                    if curr_3UTR_event_id not in All_samples_extracted_3UTR_coverage_dict_minus:
                        All_samples_extracted_3UTR_coverage_dict_minus[curr_3UTR_event_id] = []
                    All_samples_extracted_3UTR_coverage_dict_minus[curr_3UTR_event_id].append([extracted_coverage,extracted_3UTR_region])
    
    try:
        for curr_wig_file in All_Wig_files_polya:
            os.remove(tmp_dir+curr_wig_file)
    except OSError:
        pass
    else:  
        UTR_events_dict_fin = {}
        for curr_3UTR_id in UTR_events_dict:
            curr_3UTR_all_samples_bp_coverage = []
            if UTR_events_dict[curr_3UTR_id][-2] == '-' and curr_3UTR_id in All_samples_extracted_3UTR_coverage_dict_minus:
                curr_3UTR_coverage_wig = All_samples_extracted_3UTR_coverage_dict_minus[curr_3UTR_id]
                for curr_sample_curr_3UTR_coverage_wig in curr_3UTR_coverage_wig: 
                    bp_coverage = np.zeros(curr_sample_curr_3UTR_coverage_wig[1][-1] - curr_sample_curr_3UTR_coverage_wig[1][0])
                    relative_start = curr_sample_curr_3UTR_coverage_wig[1][0]
                    for i in range(len(curr_sample_curr_3UTR_coverage_wig[0])):
                        curr_region_start = curr_sample_curr_3UTR_coverage_wig[1][i] - relative_start
                        curr_region_end = curr_sample_curr_3UTR_coverage_wig[1][i+1] - relative_start
                        bp_coverage[curr_region_start:curr_region_end] = curr_sample_curr_3UTR_coverage_wig[0][i]
                    bp_coverage = bp_coverage[::-1]
                    curr_3UTR_all_samples_bp_coverage.append(bp_coverage)
            elif UTR_events_dict[curr_3UTR_id][-2] == '+' and curr_3UTR_id in All_samples_extracted_3UTR_coverage_dict:
                curr_3UTR_coverage_wig = All_samples_extracted_3UTR_coverage_dict[curr_3UTR_id]
                for curr_sample_curr_3UTR_coverage_wig in curr_3UTR_coverage_wig: 
                    bp_coverage = np.zeros(curr_sample_curr_3UTR_coverage_wig[1][-1] - curr_sample_curr_3UTR_coverage_wig[1][0])
                    relative_start = curr_sample_curr_3UTR_coverage_wig[1][0]
                    for i in range(len(curr_sample_curr_3UTR_coverage_wig[0])):
                        curr_region_start = curr_sample_curr_3UTR_coverage_wig[1][i] - relative_start
                        curr_region_end = curr_sample_curr_3UTR_coverage_wig[1][i+1] - relative_start
                        bp_coverage[curr_region_start:curr_region_end] = curr_sample_curr_3UTR_coverage_wig[0][i]
                    curr_3UTR_all_samples_bp_coverage.append(bp_coverage)
            UTR_events_dict_fin[curr_3UTR_id] = curr_3UTR_all_samples_bp_coverage
        
        min_UTR_events_dict_fin = {}
        for curr_event in UTR_events_dict_fin:
            for i in range(len(UTR_events_dict_fin[curr_event])):
                indicator_list = []
                cov_array = UTR_events_dict_fin[curr_event][i]
                probe_region = cov_array[upstream_region_one+upstream_region_two:len(cov_array)-(downstream_region_one+downstream_region_two)]
                probe_region[probe_region < min_cov] = 0
                probe_region[probe_region >= min_cov] = 1
                indicator_list.append(probe_region)
            min_UTR_events_dict_fin[curr_event] = indicator_list
            
    return min_UTR_events_dict_fin


labels = load_polya_seq_coverage(All_Wig_files_polya, 'UTR_duplicate_subtracted.bed', temp_dir)


def load_rna_seq_coverage(All_Wig_files, UTR_Annotation_file, tmp_dir):
    
    UTR_events_dict = {}
    All_Samples_Total_depth = []
    All_samples_extracted_3UTR_coverage_dict = {}
    All_samples_extracted_3UTR_coverage_dict_minus = {}
    upstream_region_one = 70
    upstream_region_two = 70
    downstream_region_one = 30
    downstream_region_two = 30
    region_chunks = 5
    
    with open(tmp_dir+UTR_Annotation_file,'r') as f: 
        for line in f:
            fields = line.strip('\n').split('\t')
            curr_chr = fields[0]
            if not curr_chr.startswith('chr'):
                curr_chr = 'chr' + curr_chr
            region_start = fields[1]
            region_end = fields[2]
            UTR_pos = '%s:%s-%s' % (curr_chr, region_start, region_end)
            UTR_events_dict[fields[3]] = [fields[0], region_start, region_end, fields[5], UTR_pos]
       
    for curr_wig_file in All_Wig_files:
        
        cur_sample_total_depth = 0
        num_line = 0
        curr_sample_All_chroms_coverage_dict = {}
        curr_sample_All_chroms_coverage_dict_minus = {}
        
        with open(tmp_dir+curr_wig_file, 'r') as f:
            for line in f:
                if '#' not in line:
                    fields = line.strip('\n').split('\t')
                    chrom_name = fields[0]
                    if not chrom_name.startswith('chr'):
                          chrom_name = 'chr' + chrom_name
                    region_start = int(fields[1])
                    region_end = int(fields[2])
                    cur_sample_total_depth += int(float(fields[-2])) * (region_end - region_start) + int(float(fields[-1])) * (region_end - region_start)
                    
                    if chrom_name not in curr_sample_All_chroms_coverage_dict:
                        curr_sample_All_chroms_coverage_dict[chrom_name] = [[0],[0]]
                    if region_start > curr_sample_All_chroms_coverage_dict[chrom_name][0][-1]:
                        curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_start)
                        curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)
                    curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_end)
                    curr_sample_All_chroms_coverage_dict[chrom_name][1].append(int(fields[-2]))
            
                    if chrom_name not in curr_sample_All_chroms_coverage_dict_minus:
                        curr_sample_All_chroms_coverage_dict_minus[chrom_name] = [[0],[0]]
                    if region_start > curr_sample_All_chroms_coverage_dict_minus[chrom_name][0][-1]:
                        curr_sample_All_chroms_coverage_dict_minus[chrom_name][0].append(region_start)
                        curr_sample_All_chroms_coverage_dict_minus[chrom_name][1].append(0)
                    curr_sample_All_chroms_coverage_dict_minus[chrom_name][0].append(region_end)
                    curr_sample_All_chroms_coverage_dict_minus[chrom_name][1].append(int(fields[-1]))
                num_line += 1
        All_Samples_Total_depth.append(cur_sample_total_depth)
        curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)
        curr_sample_All_chroms_coverage_dict_minus[chrom_name][1].append(0)
        
        for curr_3UTR_event_id in UTR_events_dict:
            curr_3UTR_structure = UTR_events_dict[curr_3UTR_event_id]
            curr_chr = curr_3UTR_structure[0]
            
            if UTR_events_dict[curr_3UTR_event_id][-2] == '+':
                if curr_chr in curr_sample_All_chroms_coverage_dict:
                    curr_chr_coverage = curr_sample_All_chroms_coverage_dict[curr_chr]
                    region_start = int(curr_3UTR_structure[1])
                    region_end = int(curr_3UTR_structure[2])
                    left_region_index = bisect(curr_chr_coverage[0],region_start)
                    right_region_index = bisect(curr_chr_coverage[0],region_end)
                    extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]
                    extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                    extracted_3UTR_region.insert(0,region_start)
                    extracted_3UTR_region.append(region_end)
                    
                    if curr_3UTR_event_id not in All_samples_extracted_3UTR_coverage_dict:
                        All_samples_extracted_3UTR_coverage_dict[curr_3UTR_event_id] = []
                    All_samples_extracted_3UTR_coverage_dict[curr_3UTR_event_id].append([extracted_coverage,extracted_3UTR_region])
           
            elif UTR_events_dict[curr_3UTR_event_id][-2] == '-':
                if curr_chr in curr_sample_All_chroms_coverage_dict_minus:
                    curr_chr_coverage = curr_sample_All_chroms_coverage_dict_minus[curr_chr]
                    region_start = int(curr_3UTR_structure[1])
                    region_end = int(curr_3UTR_structure[2])
                    left_region_index = bisect(curr_chr_coverage[0],region_start)
                    right_region_index = bisect(curr_chr_coverage[0],region_end)
                    extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]
                    extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                    extracted_3UTR_region.insert(0,region_start)
                    extracted_3UTR_region.append(region_end)
                    
                    if curr_3UTR_event_id not in All_samples_extracted_3UTR_coverage_dict_minus:
                        All_samples_extracted_3UTR_coverage_dict_minus[curr_3UTR_event_id] = []
                    All_samples_extracted_3UTR_coverage_dict_minus[curr_3UTR_event_id].append([extracted_coverage,extracted_3UTR_region])
    
    try:
        os.remove(tmp_dir+UTR_Annotation_file)
        for curr_wig_file in All_Wig_files:
            os.remove(tmp_dir+curr_wig_file)
    except OSError:
        pass
    else:  
        UTR_events_dict_fin = {}
        for curr_3UTR_id in UTR_events_dict:
            curr_3UTR_all_samples_bp_coverage = []
            if UTR_events_dict[curr_3UTR_id][-2] == '-' and curr_3UTR_id in All_samples_extracted_3UTR_coverage_dict_minus:
                curr_3UTR_coverage_wig = All_samples_extracted_3UTR_coverage_dict_minus[curr_3UTR_id]
                for curr_sample_curr_3UTR_coverage_wig in curr_3UTR_coverage_wig: 
                    bp_coverage = np.zeros(curr_sample_curr_3UTR_coverage_wig[1][-1] - curr_sample_curr_3UTR_coverage_wig[1][0])
                    relative_start = curr_sample_curr_3UTR_coverage_wig[1][0]
                    for i in range(len(curr_sample_curr_3UTR_coverage_wig[0])):
                        curr_region_start = curr_sample_curr_3UTR_coverage_wig[1][i] - relative_start
                        curr_region_end = curr_sample_curr_3UTR_coverage_wig[1][i+1] - relative_start
                        bp_coverage[curr_region_start:curr_region_end] = curr_sample_curr_3UTR_coverage_wig[0][i]
                    bp_coverage = bp_coverage[::-1]
                    curr_3UTR_all_samples_bp_coverage.append(bp_coverage)
            elif UTR_events_dict[curr_3UTR_id][-2] == '+' and curr_3UTR_id in All_samples_extracted_3UTR_coverage_dict:
                curr_3UTR_coverage_wig = All_samples_extracted_3UTR_coverage_dict[curr_3UTR_id]
                for curr_sample_curr_3UTR_coverage_wig in curr_3UTR_coverage_wig: 
                    bp_coverage = np.zeros(curr_sample_curr_3UTR_coverage_wig[1][-1] - curr_sample_curr_3UTR_coverage_wig[1][0])
                    relative_start = curr_sample_curr_3UTR_coverage_wig[1][0]
                    for i in range(len(curr_sample_curr_3UTR_coverage_wig[0])):
                        curr_region_start = curr_sample_curr_3UTR_coverage_wig[1][i] - relative_start
                        curr_region_end = curr_sample_curr_3UTR_coverage_wig[1][i+1] - relative_start
                        bp_coverage[curr_region_start:curr_region_end] = curr_sample_curr_3UTR_coverage_wig[0][i]
                    curr_3UTR_all_samples_bp_coverage.append(bp_coverage)
            UTR_events_dict_fin[curr_3UTR_id] = curr_3UTR_all_samples_bp_coverage
    
    norm_UTR_events_dict_fin = {}
    for curr_event in UTR_events_dict_fin:
        for i in range(len(UTR_events_dict_fin[curr_event])):
            exp_list = []
            exper_wide_norm = []
            cov_array = UTR_events_dict_fin[curr_event][i]
            probe_region = cov_array[upstream_region_one+upstream_region_two:len(cov_array)-(downstream_region_one+downstream_region_two)]
            exper_wide_norm.append((sum(probe_region)/np.array(All_Samples_Total_depth)[i])/len(probe_region))
            for base in range(len(cov_array)):
                 if base >= upstream_region_one+upstream_region_two and base < (len(cov_array) - (downstream_region_one+downstream_region_two)): 
                     reg_one_exp = np.sum(cov_array[base-(upstream_region_one+upstream_region_two):base-upstream_region_two])/sum(probe_region)
                     up_reg_exp_probe = [int(upstream_region_two/region_chunks)] * region_chunks
                     error = int(upstream_region_two/region_chunks * region_chunks - sum(up_reg_exp_probe))
                     up_reg_exp_probe[-1] = up_reg_exp_probe[-1] + error
                     small_up_reg = []
                     for z in range(len(up_reg_exp_probe)):
                         if z != 5:
                             small_up_reg.append(np.sum(cov_array[base-upstream_region_one+sum(up_reg_exp_probe[:z]):base-upstream_region_one+sum(up_reg_exp_probe[:z+1])])/sum(probe_region))
                     reg_two_exp = np.sum(cov_array[base+downstream_region_one:base+downstream_region_one+downstream_region_two])/sum(probe_region)
                     down_reg_exp_probe = [int(downstream_region_one/region_chunks)] * region_chunks
                     error = int(downstream_region_one/region_chunks * region_chunks - sum(down_reg_exp_probe))
                     down_reg_exp_probe[-1] = down_reg_exp_probe[-1] + error
                     small_down_reg = []
                     for z in range(len(down_reg_exp_probe)+1):
                         if z != 0:
                             small_down_reg.append(np.sum(cov_array[base+sum(down_reg_exp_probe[:z-1]):base+sum(down_reg_exp_probe[:z])])/sum(probe_region))
                     exp_list.append([np.array(reg_one_exp), np.array(small_up_reg), np.array(small_down_reg), np.array(reg_two_exp)])
            norm_UTR_events_dict_fin[curr_event] = exp_list
            norm_UTR_events_dict_fin[curr_event].append(np.array(exper_wide_norm[i]/len(cov_array)))
    
    return norm_UTR_events_dict_fin
  
                          
exp_feat = load_rna_seq_coverage(All_Wig_files, 'UTR_duplicate_subtracted.bed', temp_dir) 
    

def combine_feature_and_label_dicts(*arg):

    comb_dicts = {k:[d[k] for d in arg] for k in arg[0]}
    
    region_chunks = 5
    arr_dim = region_chunks*2+2
    fin_arr_all = np.array([]).reshape(0, 81)
    
    for event in comb_dicts:
        fin_arr = np.array([]).reshape(0, arr_dim)
        for features in range(len(comb_dicts[event])):
            if features == 0:
                for feature in range(len(comb_dicts[event][features])):
                    if feature != 1800:
                        feat_arr = np.hstack([x for x in comb_dicts[event][features][feature]])                        
                        fin_arr = np.vstack((fin_arr, feat_arr))
                    else:
                        overall_exp = np.repeat(comb_dicts[event][features][feature], len(comb_dicts[event][features])-1)
                fin_arr = np.hstack((fin_arr, overall_exp.reshape(overall_exp.shape[0], -1)))                   
            elif features == 1:
                feat_arr_fin = np.array([]).reshape(0, 64)
                for feature in range(len(comb_dicts[event][features])):
                    feat_arr = comb_dicts[event][features][feature].reshape(1, -1)
                    feat_arr_fin = np.vstack((feat_arr_fin, feat_arr))
                fin_arr = np.hstack((fin_arr, feat_arr_fin))
            elif features == 2 or features == 3 or features == 4:
                feat_arr_fin = np.array([]).reshape(0, 1)
                for feature in range(len(comb_dicts[event][features])):
                    feat_arr = np.array(int(comb_dicts[event][features][feature]))
                    feat_arr_fin = np.vstack((feat_arr_fin, feat_arr))
                fin_arr = np.hstack((fin_arr, feat_arr_fin))
            else:
                for label in range(len(comb_dicts[event][features])):
                    lab_arr = np.vstack([x for x in comb_dicts[event][features][label]])
                fin_arr = np.hstack((fin_arr, lab_arr))
        fin_arr_all = np.vstack((fin_arr_all, fin_arr))
        
    return fin_arr_all


features_and_labels = combine_feature_and_label_dicts(exp_feat, seq_feat, ssv, wsv, ca, labels)


def ml(feat_and_lab_arr):
    
    X_train_raw, X_test_raw, y_train_raw, y_test = train_test_split(feat_and_lab_arr[:, 0:80], feat_and_lab_arr[:, -1].astype('int32'), random_state = 0)
    
    scaler = MinMaxScaler()
    X_train_raw_scaled = scaler.fit_transform(X_train_raw)
    X_test = scaler.transform(X_test_raw)
    
    nm3 = NearMiss(version = 3, random_state = 0)
    X_train_unshuffled, y_train_unshuffled = nm3.fit_sample(X_train_raw_scaled, y_train_raw)
    X_train, empty_1, y_train, empty_2 = train_test_split(X_train_unshuffled, y_train_unshuffled, random_state = 0, test_size = 0)
    
    c_range = [0.1, 0.5, 1, 10, 100]
    max_iter_range = [10000, 100000]
    solver_options = ['liblinear', 'saga']
    penalty_options = ['l1', 'l2']
    param_grid = dict(C = c_range, max_iter = max_iter_range, solver = solver_options, penalty = penalty_options)
    
    f_two_scorer = make_scorer(fbeta_score, beta = 2)
    f_half_scorer = make_scorer(fbeta_score, beta = 0.5)
    scoring = ['f1', f_two_scorer, f_half_scorer, 'roc_auc', 'average_precision']
    
    lr = LogisticRegression(random_state = 0)
    
    fin_dict = {}
    for score in scoring:
        lr_grid = GridSearchCV(lr, param_grid = param_grid, scoring = score, cv = 5)
        lr_grid.fit(X_train, y_train)
        fin_dict['%s_optimized_nm3_model' % (score)] = lr_grid
    fin_dict['downsampled_nearmiss_data'] = ([X_train, X_test, y_train, y_test])
    
    return fin_dict


ml_model_and_data = ml(features_and_labels)
file = open(temp_dir+'ml_results.p', 'wb')
pickle.dump(ml_model_and_data, file)
file.close()
