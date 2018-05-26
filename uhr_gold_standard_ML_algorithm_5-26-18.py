#!/bin/bash python3

import re
import os
import datetime
import shutil
from bisect import bisect

temp_dir = '/Users/ryan/Desktop/ml_prototype/'
gold_standard_rna_seq_sorted_bam_file = temp_dir+'sorted_uhr_exp2_auto_ensembl_v91_genome_index_default_hisat2.bam'
gold_standard_polya_seq_sorted_bam_file = temp_dir+'sorted_all_uhr_polya_trimmed_ensembl_v91_genome_index_default_hisat2.bam'
gold_standard_gtf_file_from_rna_seq_data = temp_dir+'uhr_exp2_auto_ensembl_v91_genome_index_default_GRCh38.v91.txt'
chrom_sizes_file = temp_dir+'hg38_self_built.chrom.sizes'
fasta_file = temp_dir+'genome_edited_headers.fa'
sequence_bed_file = 'modified_bed_with_sequence.bed'
UTR_bed_file = 'UTR_seq.bed'

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
    
    output_write = open(tmp_dir+'UTR_seq.bed', 'w')
    
    with open(tmp_dir+sequence_bed, 'r') as f:
        for line in f:
            fields = line.strip('\n').split('\t')
            strand = fields[5]
            if strand == '+' and int(fields[10].strip(',').split(',')[-1]) >= 10000:
                all_seq = fields[-1]
                seq = all_seq[-10000:]
                tx_end = str(int(fields[1]) + 10000)
                fields[10] = fields[10].strip(',').split(',')[-1]
                fields[11] = fields[11].strip(',').split(',')[-1]
            elif strand == '+' and int(fields[10].strip(',').split(',')[-1]) < 10000:
                seq = fields[-1]    
                tx_end = str(int(fields[1]) + len(seq))
                fields[10] = fields[10].strip(',').split(',')[-1]
                fields[11] = fields[11].strip(',').split(',')[-1]
            elif strand == '-' and int(fields[10].strip(',').split(',')[0]) >= 10000:
                all_seq = fields[-1]
                seq = all_seq[-10000:]
                tx_end = str(int(fields[1]) + 10000)
                fields[10] = fields[10].strip(',').split(',')[0]
                fields[11] = fields[11].strip(',').split(',')[0]
            elif strand == '-' and int(fields[10].strip(',').split(',')[0]) < 10000:
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
    
    with open(temp_dir+UTR_bed_path, 'r') as f:
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
    return UTR_dict


modified_subtracted_utr_dict = dinuc_seq_freq('UTR_duplicate_subtracted.bed', temp_dir)
    