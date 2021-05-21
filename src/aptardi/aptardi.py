import re
import os
import sys
import datetime
from bisect import bisect
import numpy as np
import pickle
import shutil
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import LSTM, Dense, TimeDistributed, Masking, Bidirectional
from tensorflow.keras.preprocessing.sequence import pad_sequences
from tensorflow.keras.metrics import Precision, Recall
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.utils import class_weight
from tensorflow.keras.callbacks import ModelCheckpoint
import pandas as pd
import argparse
import warnings

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
warnings.filterwarnings('ignore')

def parse_args(): 
    def dir_path(outdir):
        if outdir[-1] != '/':
            outdir = outdir + '/'
        tmp_dir = datetime.datetime.now().strftime('%y%m%d_%H%M%S')
        if tmp_dir[-1] != '/':
            tmp_dir = tmp_dir + '/'
        tmp_dir = outdir+tmp_dir
        return outdir, tmp_dir
    def check_prob(value):
        if value == 'default':
            return value
        else:
            try: 
                ivalue = float(value.split('/')[0])/float(value.split('/')[1])
            except:
                try:
                    ivalue = float(value)
                except:
                    raise argparse.ArgumentTypeError('%s cannot be converted to a float' % value)
            if ivalue <= 0 or ivalue >= 1:
                raise argparse.ArgumentTypeError('%s not in (0, 1) bound' % value)
            return ivalue 
    def check_file(in_path):
        if not os.path.isfile(in_path):
            raise argparse.ArgumentTypeError('%s does not exist' % in_path)
        return in_path
    def check_positive(value):
        if value == 'default':
            return value
        else:
            try:
                ivalue = int(float(value))
                if ivalue < 0:
                    raise argparse.ArgumentTypeError('%s is not a positive integer' % value)
            except:
                 raise argparse.ArgumentTypeError('%s is not a positive integer' % value)
        return ivalue
    parser = argparse.ArgumentParser(description = 'aptardi')
    parser.add_argument('--version', '-v', action = 'version', version = 'aptardi %s' %(soft_ver))
    parser.add_argument('--o', type = dir_path, help = 'output directory', required = True, metavar = 'path')
    parser.add_argument('--r', type = argparse.FileType('r'), help = 'transcriptome reconstruction', required = True, default = sys.stdin, metavar = 'gtf/gff file')
    parser.add_argument('--f', type = check_file, help = 'fasta where headers are chromosomes', required = True, metavar = 'fasta file')
    parser.add_argument('--b', type = check_file, help = 'sorted bam file of aligned RNA-Seq reads', required = True, metavar = 'sorted bam file')
    parser.add_argument('--a', '-a', type = str, help = 'upstream/downstream mate orientations for paired-end alignment against the forward reference strand', choices = ['fr', 'rf'], required = False, default = 'fr', const = 'fr', nargs = '?', metavar = 'either fr or rf')
    parser.add_argument('--g', '-g', type = str, help = 'output gtf file', required = False, default = sys.stdout, metavar = 'output gtf')
    parser.add_argument('--d', '-d', default = False, action = 'store_true', required = False, help = 'turn debugging mode on')
    parser.add_argument('--verbose', '-vb', default = False, action = 'store_true', required = False, help = 'turn verbose mode on')
    parser.add_argument('--m', '-m', default = False, action = 'store_true', required = False, help = 'turn machine learning mode on')
    parser.add_argument('--s', '-s', type = str, help = "tab separated polyA sites file (if building model)", required = False, nargs = '?', metavar = 'polyA sites file')
    parser.add_argument('--e', '-e', type = str, help = "name to save model (if building model)", required = False, nargs = '?', metavar = 'output model') 
    parser.add_argument('--k', '-k', type = str, help = "name to save scale (if building model)", required = False, nargs = '?', metavar = 'output scale')
    parser.add_argument('--l', '-l', type = str, help = "0-based coordinates of chromosome, strand, and polyA site columns in polyA sites file (if building model)", required = False, nargs = '?', metavar = 'chromosome,strand,site')
    parser.add_argument('--p', '-p', type = check_prob, help = "probability threshold, predictions >= value will be classified as polyA sites (otherwise default used)", default = 'default', required = False, const = 'default', nargs = '?', metavar = "value in (0, 1)")
    parser.add_argument('--c', '-c', type = int, help = "seed state to set for train test split if building model (default no seed)", required = False, nargs = '?', metavar = 'integer')
    parser.add_argument('--n', '-n', type = str, help = "pre-built model (if not building model)", required = False, nargs = '?', metavar = 'input model')
    parser.add_argument('--t', '-t', type = str, help = "scale from pre-built model (if not building model)", required = False, nargs = '?', metavar = 'input scale')
    parser.add_argument('--i', '-i', type = check_positive, help = "maximum number of bins analyzed per transcript (if building model)", required = False, default = 300, const = 300, nargs = '?', metavar = 'positive integer (default 300) total max length of transcript analyzed is this * bin size')
    parser.add_argument('--w', '-w', type = int, help = 'size of bins for making predictions', choices = [25, 50, 75, 100, 125, 150, 175, 200], required = False, default = 100, const = 100, nargs = '?', metavar = 'choose between 25-200 in 25 base increments (default 100)')
    return parser.parse_args()

def check_ml_file(arg_in, tmp_dir, outdir):
    if machine_learning_mode:
        if arg_in.n or arg_in.t:
            sys.stderr.write('If machine learning mode, please do not provide any of the following arguments: pre-built model (-n) or pre-built model scale (-t)'+'\n')
            shutil.rmtree(tmp_dir, ignore_errors = True)
            try:
                os.rmdir(outdir)
            except OSError:
                pass
            sys.exit(1)  
        if not arg_in.s or not os.path.isfile(arg_in.s):
            sys.stderr.write('If machine learning mode, a tab delimited file for labels is required'+'\n')
            shutil.rmtree(tmp_dir, ignore_errors = True)
            try:
                os.rmdir(outdir)
            except OSError:
                pass
            sys.exit(1)  
        if not arg_in.e or not arg_in.e.endswith('.hdf5'):
            sys.stderr.write('If machine learning mode, a name to save the model as, with .hdf5 extension, is required'+'\n')
            shutil.rmtree(tmp_dir, ignore_errors = True)
            try:
                os.rmidr(outdir)
            except OSError:
                pass
            sys.exit(1)
        if not arg_in.k or not arg_in.k.endswith('.pk'):
            sys.stderr.write('If machine learning mode, a name to save the model scale as, with .pk extension, is required'+'\n')
            shutil.rmtree(tmp_dir, ignore_errors = True)
            try:
                os.rmdir(outdir)
            except OSError:
                pass
            sys.exit(1)         
        try:
            chrom_lab, strand_lab, site_lab = [int(float(x)) for x in arg_in.l.strip().split(',')]
        except:
            sys.stderr.write('If machine learning mode, a comma separated list denoting chromosome, strand, and polyA site columns in tab delimited label file is also required'+'\n')
            shutil.rmtree(tmp_dir, ignore_errors = True)
            try:
                os.rmdir(outdir)
            except OSError:
                pass
            sys.exit(1)  
        return chrom_lab, strand_lab, site_lab, arg_in.s, arg_in.e, arg_in.k, arg_in.i
    else:
        return None, None, None, None, None, None, None
    
def check_model_files(arg_in, tmp_dir, outdir):
    if not machine_learning_mode:
        if arg_in.s or arg_in.e or arg_in.k or arg_in.l:
            sys.stderr.write('Since not machine learning mode, please do not provide any of the following arguments: label file (-s), columns of label file (-l), name to save model (-e), or name to save model scale (-k)'+'\n')
            shutil.rmtree(tmp_dir, ignore_errors = True)
            try:
                os.rmdir(outdir)
            except OSError:
                pass
            sys.exit(1)
        if not arg_in.n or not os.path.isfile(arg_in.n):
            sys.stderr.write('If not machine learning mode, pre-built model required'+'\n')
            shutil.rmtree(tmp_dir, ignore_errors = True)
            try:
                os.rmdir(outdir)
            except OSError:
                pass
            sys.exit(1)  
        if not arg_in.t or not os.path.isfile(arg_in.t):
            sys.stderr.write('If not machine learning mode, scale from pre-built model required'+'\n')
            shutil.rmtree(tmp_dir, ignore_errors = True)
            try:
                os.rmdir(outdir)
            except OSError:
                pass
            sys.exit(1) 
        return arg_in.n, arg_in.t
    else:
        return None, None
       
def time_now(proc, starting = False, finishing = False):
    curr_time = datetime.datetime.now()
    time = curr_time.strftime('%d %b %Y %I:%M %p')
    if starting:
        sys.stderr.write(time+': '+proc+'\n')
    elif finishing:
        sys.stderr.write('Finished! ('+ time+')'+'\n') 

def stg_gtf_to_bed(gtf_file, tmp_dir):
    def bed_line(estart_val, eend_val, field_vals, nline_val, fpkm_val, tpm_val, coverage_val):
        estp = estart_val[0] - 1
        eedp = eend_val[-1]
        allids = {}
        curr_att_field = field_vals[8]
        curr_strand = field_vals[6]
        curr_chrom = field_vals[0]
        geneid = re.findall(r'gene_id \"([\w\.]+)\"', curr_att_field) 
        if len(geneid) == 0:
            geneid = 'NA'
        else:
            geneid = geneid[0]
        refgeneid = re.findall(r'ref_gene_id \"([\w\.]+)\"', curr_att_field)
        if len(refgeneid) == 0:
            refgeneid = 'NA'
        else:
            refgeneid = refgeneid[0]
        refgenename = re.findall(r'ref_gene_name \"([\w\.]+)\"', curr_att_field)
        if len(refgenename) == 0:
            refgenename = 'NA'
        else:
            refgenename = refgenename[0]  
        transid = re.findall(r'transcript_id \"([\w\.]+)\"', curr_att_field) 
        if len(transid) == 0:
            transid = 'Trans_' + str(nline_val)
        else:
            transid = transid[0]
        if transid in allids.keys():
            transid2 = transid + '_DUP' + str(allids[transid])
            allids[transid] = allids[transid] + 1
            transid = transid2
        else:
            allids[transid] = 1  
        reftranscriptid = re.findall(r'reference_id \"([\w\.]+)\"', curr_att_field) 
        if len(reftranscriptid) == 0:
            reftranscriptid = 'NA'
        else:
            reftranscriptid = reftranscriptid[0]
        if not curr_chrom.startswith('chr'):
            curr_chrom = 'chr' + curr_chrom
        try:
            if curr_strand == '+':
                cov_fin = float(coverage_val[-1])
            elif curr_strand == '-':
                cov_fin = float(coverage_val[0])
        except:
            cov_fin = 'NA'
        line = [curr_chrom, str(estp), str(eedp), '|'.join([transid, geneid, curr_chrom, str(estp), str(eedp), curr_strand]), '|'.join([reftranscriptid, refgeneid, refgenename]), curr_strand, fpkm_val, tpm_val, str(cov_fin), str(len(estart_val))]
        seglen = [eend_val[i] - estart_val[i] + 1 for i in range(len(estart_val))]
        segstart = [estart_val[i] - estart_val[0] for i in range(len(estart_val))]
        strl = str(seglen[0])
        for i in range(1, len(seglen)):
            strl += ',' + str(seglen[i])
        strs = str(segstart[0])
        for i in range(1, len(segstart)):
            strs +=',' + str(segstart[i])
        line.extend([strl, strs])
        return line
    if verbose:
        time_now('Loading transcript reconstructions..', starting = True)  
    extracted_bed_line = []   
    estart = []
    eend = []
    fpkm = ''
    tpm = ''
    cov = []
    nline = 0
    prevfield = []
    prevtransid = ''
    for lines in gtf_file:
        if lines.strip().split('\t')[0] == '\n' or lines[0] == '#' or not lines.strip():
            continue
        fields = lines.strip().split('\t')
        fields[8:] = [(' ').join(fields[8:])]
        att_field = fields[8]
        type_field = fields[2]
        recon_field = fields[1]
        chrom = fields[0]
        strand = fields[6]
        if '_' in chrom or strand not in ['-', '+']:
            continue
        nline = nline + 1
        if nline == 1 and recon_field not in ['Cufflinks', 'StringTie', 'aptardi'] and verbose:
            sys.stderr.write('Warning: The GTF/GFF reconstruction file may not have been produced by Cufflinks or StringTie or aptardi and therefore may give erroneous results'+'\n')
        if type_field != 'exon' and type_field != 'transcript':
            continue        
        transid = [att_field.split('transcript_id ')[1].split()[0].replace(';', '').replace('"', '')]
        if len(transid) > 0:
            transid = transid[0]
        else:
            transid = ''
        if type_field == 'transcript' or (prevtransid != '' and transid != '' and transid != prevtransid):
            if len(estart) != 0:
                order_line = 0
                for a, b in zip(estart, estart[1:]):
                    if a > b:
                        order_line += 1
                if order_line == len(estart)-1 and len(estart) > 1:
                    extracted_bed_line.append(bed_line(estart[::-1], eend[::-1], prevfield, nline, fpkm, tpm, cov[::-1]))
                else:
                   extracted_bed_line.append(bed_line(estart, eend, prevfield, nline, fpkm, tpm, cov))
            estart = []
            eend = []
            cov = []
            fpkm = ''
            tpm = ''
        prevfield = fields
        prevtransid = transid
        if type_field == 'exon':
            est = int(float(fields[3]))
            eed = int(float(fields[4]))
            try:
                covval = float(re.findall(r'cov \"([\w\.]+)\"', att_field)[0])
            except:
                covval = 'NA'
            estart += [est]
            eend += [eed]
            cov += [covval]
        if type_field == 'transcript':
            fpkmval = re.findall(r'FPKM \"([\w\.]+)\"', att_field)
            tpmval = re.findall(r'TPM \"([\w\.]+)\"', att_field)
        if len(fpkmval) == 0:
            fpkm = '0'
        else:
            fpkm = fpkmval[0]
        if len(tpmval) == 0:
            tpm = '0'
        else:
            tpm = tpmval[0]
    if len(estart) != 0 and '_' not in chrom and strand in ['-', '+']:
        order_line = 0
        for a, b in zip(estart, estart[1:]):
            if a > b:
                order_line += 1
        if order_line == len(estart)-1 and len(estart) > 1:
            extracted_bed_line.append(bed_line(estart[::-1], eend[::-1], prevfield, nline, fpkm, tpm, cov[::-1]))
        else:
           extracted_bed_line.append(bed_line(estart, eend, prevfield, nline, fpkm, tpm, cov))
    if len(extracted_bed_line) == 0:
        sys.stderr.write('Error: No transcript reconstructions identified, perhaps the GTF/GFF reconstruction file format does not match that produced by StringTie/Cufflinks?'+'\n')
        if not debugging:
            shutil.rmtree(tmp_dir, ignore_errors = True)
            try:
                os.rmdir(out_dir)
            except OSError:
                pass
        sys.exit(1)
    return extracted_bed_line

def format_bed_list(b_list):
    fin_list = []
    for i in b_list:
        if (not str.isdigit(i[3].split('|')[2][3:]) and not i[3].split('|')[2][3:].upper() == 'X' and not i[3].split('|')[2][3:].upper() == 'Y') or int(float(i[3].split('|')[3])) < 0:
            continue
        start = int(float(i[1]))
        size = [int(float(x)) for x in i[10].split(',')]
        pos = [int(float(x)) for x in i[11].split(',')]
        abs_pos_start = [x + start for x in pos]
        abs_pos_end = [sum(x) for x in zip(abs_pos_start, size)]
        junc = [val for pair in zip(abs_pos_start, abs_pos_end) for val in pair]
        fin_list.append([i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10], i[11], (',').join([str(a) for a in junc])])
    if verbose:
        time_now('', finishing = True) 
    return fin_list

def bam_to_bedgraph(rna_seq_bam_path, tmp_dir, save_file):
    if verbose:
        time_now('Converting bam file to bedgraph file..', starting = True)   
    os.system('samtools view -h -f 0X40 %s | bedtools genomecov -bg -split -strand + -ibam - > %s' % (rna_seq_bam_path, tmp_dir+'first_plus.bg'))
    os.system('samtools view -h -f 0X40 %s | bedtools genomecov -bg -split -strand - -ibam - > %s' % (rna_seq_bam_path, tmp_dir+'first_minus.bg'))
    os.system('samtools view -h -f 0X80 %s | bedtools genomecov -bg -split -strand + -ibam - > %s' % (rna_seq_bam_path, tmp_dir+'second_plus.bg'))
    os.system('samtools view -h -f 0X80 %s | bedtools genomecov -bg -split -strand - -ibam - > %s' % (rna_seq_bam_path, tmp_dir+'second_minus.bg'))
    if orientation == 'fr':
        os.system('''bedtools unionbedg -i %s %s | awk '{print $1 "\t" $2 "\t" $3 "\t" ($4+$5)}' - > %s''' % (tmp_dir+'first_minus.bg', tmp_dir+'second_plus.bg', tmp_dir+'comb_plus.bg'))
        os.system('''bedtools unionbedg -i %s %s | awk '{print $1 "\t" $2 "\t" $3 "\t" ($4+$5)}' - > %s''' % (tmp_dir+'first_plus.bg', tmp_dir+'second_minus.bg', tmp_dir+'comb_minus.bg'))
        os.system('bedtools unionbedg -i %s %s > %s' % (tmp_dir+'comb_plus.bg', tmp_dir+'comb_minus.bg', tmp_dir+save_file))
    elif orientation == 'rf':
        os.system('''bedtools unionbedg -i %s %s | awk '{print $1 "\t" $2 "\t" $3 "\t" ($4+$5)}' - > %s''' % (tmp_dir+'first_minus.bg', tmp_dir+'second_plus.bg', tmp_dir+'comb_minus.bg'))
        os.system('''bedtools unionbedg -i %s %s | awk '{print $1 "\t" $2 "\t" $3 "\t" ($4+$5)}' - > %s''' % (tmp_dir+'first_plus.bg', tmp_dir+'second_minus.bg', tmp_dir+'comb_plus.bg'))
        os.system('bedtools unionbedg -i %s %s > %s' % (tmp_dir+'comb_plus.bg', tmp_dir+'comb_minus.bg', tmp_dir+save_file))
    if not debugging:
        try:
            os.remove(tmp_dir+'comb_plus.bg')
        except OSError:
            pass
        try:
            os.remove(tmp_dir+'comb_minus.bg')
        except OSError:
            pass
        try:
            os.remove(tmp_dir+'first_plus.bg')
        except OSError:
            pass
        try:
            os.remove(tmp_dir+'first_minus.bg')
        except OSError:
            pass
        try:
            os.remove(tmp_dir+'second_plus.bg')
        except OSError:
            pass
        try:
            os.remove(tmp_dir+'second_minus.bg')
        except OSError:
            pass
    if verbose:
        time_now('', finishing = True) 

def get_chrom_sizes(fa_path, chr_file, tmp_dir):
    os.system('samtools faidx %s' % (fa_path))
    os.system('cut -f1,2 %s > %s' % (fa_path+'.fai', tmp_dir+chr_file))
    os.system('rm %s' % (fa_path+'.fai'))
        
def load_chrom_sizes(chrom_file, tmp_dir):
    if verbose:
        time_now('Processing reconstructions and defining data instances..', starting = True)  
    chrom_sizes_dict = {}
    try:
        with open(tmp_dir+chrom_file, 'rt') as f:
            for line in f:
                fields = line.strip('\n').split('\t')
                if line.strip().split('\t')[0] == '\n' or line[0] == '#' or not line.strip():
                    continue
                chrom = fields[0]
                chrom_size = fields[1]
                if not chrom.startswith('chr'):
                    chrom = 'chr' + chrom
                chrom_sizes_dict[chrom] = [chrom_size]
    except:
        sys.stderr.write('Could not open chromosome sizes file'+'\n')
        if not debugging:
            shutil.rmtree(tmp_dir, ignore_errors = True)
            try:
                os.rmdir(out_dir)
            except OSError:
                pass
        sys.exit(1) 
    if not debugging:
        try:
            os.remove(tmp_dir+chrom_file)
        except OSError:
            pass
    return chrom_sizes_dict

def extract_extended_and_orig_end_utr(bed_list_raw, chrom_sizes_dict, tmp_dir, extended_end_utr_file):
    utr_size = utr_extension + window_size * 2
    output_write = open(tmp_dir+extended_end_utr_file, 'w')
    for fields in bed_list_raw:
        strand = fields[5]
        tx_start = fields[1]
        tx_end = fields[2]  
        abs_exon_starts = [int(float(tx_start))+int(float(i)) for i in fields[11].strip(',').split(',')]
        exon_sizes = [int(float(i)) for i in fields[10].strip(',').split(',')]
        if strand == '+' and int(float(fields[10].strip(',').split(',')[-1])) < utr_size:
            extended_UTR = str(int(float(tx_end)) + (utr_size - int(float(fields[10].strip(',').split(',')[-1])))) 
            if fields[0] in chrom_sizes_dict:
                total_chrom_length = chrom_sizes_dict[fields[0]][0]
                if int(float(total_chrom_length)) < int(float(extended_UTR)):
                    extended_UTR = total_chrom_length
            tx_end = extended_UTR
            exon_sizes[-1] = int(float(tx_end))-abs_exon_starts[-1]
        elif strand == '-' and int(float(fields[10].strip(',').split(',')[0])) < utr_size:
            extended_UTR = str(int(float(tx_start)) - (utr_size - int(float(fields[10].strip(',').split(',')[0]))))
            if int(float(extended_UTR)) < 0:
                extended_UTR = '0'
            tx_start = extended_UTR
            abs_exon_starts[0] = int(float(extended_UTR))
            exon_sizes[0] = int(float(fields[1]))-int(float(tx_start))+exon_sizes[0]
        elif strand == '+' and int(float(fields[10].strip(',').split(',')[-1])) >= utr_size or strand == '-' and int(float(fields[10].strip(',').split(',')[0])) >= utr_size:
            pass
        else:
            continue          
        rel_exon_starts = (',').join([str(int(float(i))-int(float(tx_start))) for i in abs_exon_starts])
        exon_sizes = (',').join(str(i) for i in exon_sizes)
        if strand == '+':
            tx_start = str(int(float(tx_end)) - int(float(exon_sizes.strip(',').split(',')[-1])))
            exon_sizes = exon_sizes.strip(',').split(',')[-1]
            rel_exon_starts = (',').join([str(int(float(i))-int(float(tx_start))) for i in abs_exon_starts])[-1]
        elif strand == '-':
            tx_end = str(int(float(tx_start)) + int(float(exon_sizes.strip(',').split(',')[0])))
            exon_sizes = exon_sizes.strip(',').split(',')[0]
            rel_exon_starts = (',').join([str(int(float(i))-int(float(tx_start))) for i in abs_exon_starts])[0]    
        else:
            continue
        write_line = [fields[0], tx_start, tx_end, fields[3], fields[7], strand, fields[6], fields[7], fields[12], '1', exon_sizes, rel_exon_starts]
        strand = fields[5]
        tx_start = fields[1]
        tx_end = fields[2]  
        abs_exon_starts = [int(float(tx_start))+int(float(i)) for i in fields[11].strip(',').split(',')]
        exon_sizes = [int(float(i)) for i in fields[10].strip(',').split(',')]         
        if strand == '-':
            tx_end = str(int(float(tx_start)) + int(float(exon_sizes[0])))
        elif strand == '+':
            tx_start = str(int(float(tx_end)) - int(float(exon_sizes[-1])))
        else:
            continue
        write_line[6] = tx_start
        write_line[7] = tx_end
        output_write.writelines('\t'.join(write_line) + '\n')
    output_write.close()

def extract_extension(bed_list_raw, chrom_sizes_dict, tmp_dir, out_file):
    utr_size = utr_extension + window_size * 2
    output_write = open(tmp_dir+out_file, 'w')
    for fields in bed_list_raw:
        strand = fields[5]
        tx_start = fields[1]
        tx_end = fields[2]  
        abs_exon_starts = [int(float(tx_start))+int(float(i)) for i in fields[11].strip(',').split(',')]
        exon_sizes = [int(float(i)) for i in fields[10].strip(',').split(',')]
        if strand == '+' and int(float(fields[10].strip(',').split(',')[-1])) < utr_size:
            extended_UTR = str(int(float(tx_end)) + (utr_size - int(float(fields[10].strip(',').split(',')[-1])))) 
            if fields[0] in chrom_sizes_dict:
                total_chrom_length = chrom_sizes_dict[fields[0]][0]
                if int(float(total_chrom_length)) < int(float(extended_UTR)):
                    extended_UTR = total_chrom_length
            tx_end = extended_UTR
            exon_sizes[-1] = int(float(tx_end))-abs_exon_starts[-1]
        elif strand == '-' and int(float(fields[10].strip(',').split(',')[0])) < utr_size:
            extended_UTR = str(int(float(tx_start)) - (utr_size - int(float(fields[10].strip(',').split(',')[0]))))
            if int(float(extended_UTR)) < 0:
                extended_UTR = '0'
            tx_start = extended_UTR
            abs_exon_starts[0] = int(float(extended_UTR))
            exon_sizes[0] = int(float(fields[1]))-int(float(tx_start))+exon_sizes[0] 
        else:
            continue          
        write_line = [fields[0], tx_start, tx_end, fields[3], fields[12], strand]
        if strand == '-':
            write_line[2] = fields[1]
        elif strand == '+':
            write_line[1] = fields[2]
        output_write.writelines('\t'.join(write_line) + '\n')

def extract_start_utr(bed_list_raw, tmp_dir, out_file):
    output_write = open(tmp_dir+out_file, 'w')
    for fields in bed_list_raw:
        strand = fields[5]
        tx_start = fields[1]
        tx_end = fields[2]  
        exon_sizes = [int(float(i)) for i in fields[10].strip(',').split(',')]         
        if strand == '+':
            tx_end = str(int(float(tx_start)) + int(float(exon_sizes[0])))
            exon_sizes = str(exon_sizes[0])
        elif strand == '-':
            tx_start = str(int(float(tx_end)) - int(float(exon_sizes[-1])))
            exon_sizes = str(exon_sizes[-1]) 
        else:
            continue
        write_line = [fields[0], tx_start, tx_end, fields[3], fields[12], strand]
        output_write.writelines('\t'.join(write_line) + '\n')
    output_write.close()

def intersect_extensions_start_utrs(extension_file, start_utr_file, tmp_dir, intersected_extension_start_utr_file):
    os.system('bedtools intersect -a %s -b %s -s -wo > %s' % (tmp_dir+extension_file, tmp_dir+start_utr_file, tmp_dir+intersected_extension_start_utr_file))
    if not debugging:
        try:
            os.remove(tmp_dir+extension_file)
        except OSError:
            pass
        try:
            os.remove(tmp_dir+start_utr_file)
        except OSError:
            pass

def refine_utr_extension(tmp_dir, intersected_extension_start_utr_file, extended_end_utr_file, refined_utr_file):
    extended_utr_subtract_overlap = {}
    overlap_extension = {}
    output_write = open(tmp_dir+refined_utr_file, 'w') 
    with open(tmp_dir+extended_end_utr_file, 'r') as f:
        for line in f:
            if line.strip().split('\t')[0] == '\n' or line[0] == '#' or not line.strip():
                continue
            fields = line.strip().split('\t')
            name = fields[3]
            extended_utr_subtract_overlap[name] = fields
    with open(tmp_dir+intersected_extension_start_utr_file, 'r') as f:
        for line in f:
            if line.strip().split('\t')[0] == '\n' or line[0] == '#' or not line.strip():
                continue
            fields = line.strip().split('\t')
            name = fields[3]
            strand = fields[5]
            start = fields[7]
            end = fields[8]
            if name not in overlap_extension:
                overlap_extension[name] = [fields[0:6], strand, [start, end]]
            else:
                overlap_extension[name][2].extend([start, end])                 
    for ev in extended_utr_subtract_overlap:
        if ev in overlap_extension:
            if overlap_extension[ev][1] == '-':
                if int(float(max(overlap_extension[ev][2]))) < int(float(extended_utr_subtract_overlap[ev][2])):
                    extended_utr_subtract_overlap[ev][1] = max(overlap_extension[ev][2])
                else:
                    extended_utr_subtract_overlap[ev][1] = extended_utr_subtract_overlap[ev][6]
            elif overlap_extension[ev][1] == '+':
                if int(float(min(overlap_extension[ev][2]))) > int(float(extended_utr_subtract_overlap[ev][1])):
                    extended_utr_subtract_overlap[ev][2] = min(overlap_extension[ev][2])
                else:
                     extended_utr_subtract_overlap[ev][2] = extended_utr_subtract_overlap[ev][7]
            extended_utr_subtract_overlap[ev][10] = str(int(float(extended_utr_subtract_overlap[ev][2])) - int(float(extended_utr_subtract_overlap[ev][1]))) 
        output_write.writelines('\t'.join(extended_utr_subtract_overlap[ev]) + '\n')
    output_write.close()
    if not debugging:
        try:
            os.remove(tmp_dir+extended_end_utr_file)
        except OSError:
            pass
        try:
            os.remove(tmp_dir+intersected_extension_start_utr_file)
        except OSError:
            pass
  
def extract_rna_seq_coverage_utr(bedgraph_files, tmp_dir, refined_utrs_file):
    fl_dict = {}
    All_Samples_Total_depth = []
    cur_sample_total_depth = 0
    curr_sample_All_chroms_coverage_dict = {}
    curr_sample_All_chroms_coverage_dict_minus = {}
    with open(tmp_dir+bedgraph_files, 'r') as f:
        for line in f:
            if line.strip().split('\t')[0] == '\n' or line[0] == '#' or not line.strip():
                continue
            fields = line.strip('\n').split('\t')
            chrom_name = fields[0]
            if not chrom_name.startswith('chr'):
                  chrom_name = 'chr' + chrom_name
            region_start = int(float(fields[1]))
            region_end = int(float(fields[2]))
            cur_sample_total_depth += int(float(fields[-2])) * (region_end - region_start) + int(float(fields[-1])) * (region_end - region_start)
            if chrom_name not in curr_sample_All_chroms_coverage_dict:
                curr_sample_All_chroms_coverage_dict[chrom_name] = [[0],[0]]
            if region_start > curr_sample_All_chroms_coverage_dict[chrom_name][0][-1]:
                curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_start)
                curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)
            curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_end)
            curr_sample_All_chroms_coverage_dict[chrom_name][1].append(int(float(fields[-2])))
            if chrom_name not in curr_sample_All_chroms_coverage_dict_minus:
                curr_sample_All_chroms_coverage_dict_minus[chrom_name] = [[0],[0]]
            if region_start > curr_sample_All_chroms_coverage_dict_minus[chrom_name][0][-1]:
                curr_sample_All_chroms_coverage_dict_minus[chrom_name][0].append(region_start)
                curr_sample_All_chroms_coverage_dict_minus[chrom_name][1].append(0)
            curr_sample_All_chroms_coverage_dict_minus[chrom_name][0].append(region_end)
            curr_sample_All_chroms_coverage_dict_minus[chrom_name][1].append(int(float(fields[-1])))
    if cur_sample_total_depth != 0:
        All_Samples_Total_depth.append(cur_sample_total_depth)
    else:
        sys.stderr.write('Error: No read coverage detected'+'\n')
        if not debugging:
            shutil.rmtree(tmp_dir, ignore_errors = True)
            try:
                os.rmdir(out_dir)
            except OSError:
                pass
        sys.exit(1)
    curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)
    curr_sample_All_chroms_coverage_dict_minus[chrom_name][1].append(0)
    with open(tmp_dir+refined_utrs_file, 'r') as f:
        for line in f:
            if line.strip().split('\t')[0] == '\n' or line[0] == '#' or not line.strip():
                continue
            fields = line.strip('\n').split('\t')
            curr_chr = fields[0]
            strand = fields[5]
            ev = fields[3]
            tmp_junc = [int(float(ex)) for ex in fields[8].split(',')]
            if strand == '+':
                tmp_junc.extend([int(float(fields[1])), int(float(fields[2]))])
            elif strand == '-':
                tmp_junc = ([int(float(fields[1])), int(float(fields[2]))]) + tmp_junc
            junc = [tmp_junc[g:g+2] for g in range(0, len(tmp_junc), 2)]
            tmp_region_start, tmp_region_end = map(list, zip(*junc))
            if ev not in fl_dict:
                fl_dict[ev] = []
            tmp_reg_list = []
            for region_start, region_end in zip(tmp_region_start, tmp_region_end):
                if strand == '+' and curr_chr in curr_sample_All_chroms_coverage_dict:
                    curr_chr_coverage = curr_sample_All_chroms_coverage_dict[curr_chr]
                    left_region_index = bisect(curr_chr_coverage[0], region_start)
                    right_region_index = bisect(curr_chr_coverage[0], region_end)
                    extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]
                    extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                    extracted_3UTR_region.insert(0, region_start)
                    extracted_3UTR_region.append(region_end)
                    bp_coverage = np.zeros(extracted_3UTR_region[-1] - extracted_3UTR_region[0])
                elif strand == '-' and curr_chr in curr_sample_All_chroms_coverage_dict_minus:
                    curr_chr_coverage = curr_sample_All_chroms_coverage_dict_minus[curr_chr]
                    left_region_index = bisect(curr_chr_coverage[0], region_start)
                    right_region_index = bisect(curr_chr_coverage[0], region_end)
                    extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]
                    extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                    extracted_3UTR_region.insert(0, region_start)
                    extracted_3UTR_region.append(region_end)
                    bp_coverage = np.zeros(extracted_3UTR_region[-1] - extracted_3UTR_region[0])
                else:
                    continue
                relative_start = extracted_3UTR_region[0]
                for i in range(len(extracted_coverage)):
                    curr_region_start = extracted_3UTR_region[i] - relative_start
                    curr_region_end = extracted_3UTR_region[i+1] - relative_start
                    bp_coverage[curr_region_start:curr_region_end] = extracted_coverage[i]
                if strand == '-':
                    bp_coverage = bp_coverage[::-1]
                tmp_reg_list.append(bp_coverage)
            if strand == '-':
                tmp_reg_list = tmp_reg_list[::-1]
            fl_dict[ev].extend(tmp_reg_list)
    return fl_dict, All_Samples_Total_depth

def identify_transcript_end_base(trans_rna_cov_dict):
    fin_dict = {}
    for utr_event, trans_cov in trans_rna_cov_dict.items():
        fin_dict[utr_event] = []
        try:
            trans_cov_summary = np.mean(np.concatenate([x for x in trans_cov[0:-1]]))
            beg_utr_cov = np.mean(trans_cov[-2][0:window_size])
            extended_utr = trans_cov[-1]
        except:
            fin_dict[utr_event] = [0]
            continue
        if not beg_utr_cov >= trans_cov_summary*start_utr_cutoff:
            fin_dict[utr_event] = [0]
            continue
        if not len(extended_utr) >= window_size*min_num_regions:
            fin_dict[utr_event] = [0]
            continue
        else:
            window_cov_cutoff = beg_utr_cov*per_of_beg_utr
            for i in range(window_size, len(extended_utr)-window_size*2+1):
                prev_cov = extended_utr[i-window_size:i]
                curr_cov = extended_utr[i:i+window_size]
                next_cov = extended_utr[i+window_size:i+window_size*2]
                bases_above_cutoff = len(curr_cov[np.where(curr_cov > window_cov_cutoff)])
                if bases_above_cutoff >= len(curr_cov)*per_bases_curr_window:
                    continue
                elif not np.mean(prev_cov) < beg_utr_cov*per_of_beg_utr and np.mean(next_cov) < beg_utr_cov*per_of_beg_utr: 
                    continue
                elif not extended_utr[i] < beg_utr_cov*per_of_beg_utr:
                    continue
                else:
                    fin_dict[utr_event] = [i]
                    break
            else:
                fin_dict[utr_event] = [i]
    return fin_dict

def redefine_refined_utrs(refined_utr_file, ending_base_dict, new_refined_utr_file, tmp_dir):
    tmp_dict = {}
    output_write = open(tmp_dir+new_refined_utr_file, 'w') 
    with open(tmp_dir+refined_utr_file, 'r') as f:
        for line in f:
            if line.strip().split('\t')[0] == '\n' or line[0] == '#' or not line.strip():
                continue
            fields = line.strip().split('\t')
            tmp_dict[fields[3]] = fields
    for ev in tmp_dict:
        if ending_base_dict[ev][0] == 0:
            continue
        strand = tmp_dict[ev][5]
        start = int(float(tmp_dict[ev][1]))
        end = int(float(tmp_dict[ev][2]))
        if strand == '+':
            tmp_dict[ev][2] = str(start + int(float(ending_base_dict[ev][0])))
        elif strand == '-':
            tmp_dict[ev][1] = str(end - int(float(ending_base_dict[ev][0])))
        tmp_dict[ev] = tmp_dict[ev][0:-3]
        tmp_dict[ev].append(str(int(float(tmp_dict[ev][2]))-int(float(tmp_dict[ev][1]))))
        output_write.writelines('\t'.join(tmp_dict[ev]) + '\n')
    output_write.close() 
    if not debugging:
        try:
            os.remove(tmp_dir+refined_utr_file)
        except OSError:
            pass

def window_refined_utrs(tmp_dir, refined_utr_file, window_refined_utr_file):
    output_write = open(tmp_dir+window_refined_utr_file, 'w') 
    with open(tmp_dir+refined_utr_file, 'r') as f:
        for line in f:
            if line.strip().split('\t')[0] == '\n' or line[0] == '#'  or not line.strip():
                continue
            fields = line.strip().split('\t')
            start = int(float(fields[1]))
            end = int(float(fields[2]))
            if end-start != utr_extension:
                utr_size = str(end-start-(end-start)%window_size+window_size)
            else:
                utr_size = str(utr_extension)
            if fields[5] == '+':
                fields[2] = str(start+int(float(utr_size)))
            elif fields[5] == '-':
                fields[1] = str(end-int(float(utr_size)))
            fields.append(utr_size)
            output_write.writelines('\t'.join(fields) + '\n')
    output_write.close() 
    if not debugging:
        try:
            os.remove(tmp_dir+refined_utr_file)
        except OSError:
            pass
        
def merge_refined_utrs(tmp_dir, refined_utr_file, merged_refined_utr_file):
    if machine_learning_mode:
        os.system('sort -k1,1 -k2,2n %s | bedtools merge -i - -s -c 4,5,6,9 -o collapse,collapse,distinct,collapse -delim "/" > %s' % (tmp_dir+refined_utr_file, tmp_dir+merged_refined_utr_file))
 
def define_data_instances(tmp_dir, instances_file, chrom_sizes_dict):
    if machine_learning_mode:
        window_refined_utrs_file = 'merged.bed'
    else:
        window_refined_utrs_file = 'window.bed'
    output_write = open(tmp_dir+instances_file, 'w')
    with open(tmp_dir+window_refined_utrs_file, 'r') as f:
        for line in f:
            if line.strip().split('\t')[0] == '\n' or line[0] == '#' or not line.strip():
                continue
            fields = line.strip('\n').split('\t')
            name = fields[3]
            start = int(float(fields[1]))
            end = int(float(fields[2]))
            strand = fields[5]
            tpm = fields[4]
            comb_tpm = str(sum([float(x) for x in tpm.split('/')]))
            for i in range(start, end, window_size):
                base_start = i - dna_probe_region_size
                base_end = i + dna_probe_region_size + window_size
                if base_start < 0:
                    continue
                if fields[0] in chrom_sizes_dict:
                    total_chrom_length = chrom_sizes_dict[fields[0]][0]
                    if int(float(total_chrom_length)) < base_end:
                        continue
                write_line = [fields[0], str(base_start), str(base_end), name, comb_tpm, strand, str(i), str(i + window_size)]
                output_write.writelines('\t'.join(write_line) + '\n')
    output_write.close()
    if not debugging:
        try:
            os.remove(tmp_dir+window_refined_utrs_file)
        except OSError:
            pass
    if verbose:
        time_now('', finishing = True) 
  
def check_fasta_headers(local_fasta_file, tmp_dir, head_file):
    os.system('''awk 'sub(/^>/, "")' %s > %s''' % (local_fasta_file, tmp_dir+head_file))
    n = 0
    nline = 0
    with open(tmp_dir+head_file) as f:
        for line in f:
            nline += 1
            if line.startswith('chr'):
                n += 1
    if n == nline:
        new_fasta_file = None
        if not debugging:
            try:
                os.remove(tmp_dir+head_file) 
            except OSError:
                pass
        pass
    elif n == 0:
        new_fasta_file = tmp_dir+'tmp.fa'
        os.system('''sed "s/^>/>chr/" %s > %s''' % (local_fasta_file, new_fasta_file))
        if not debugging:
            try:
                os.remove(tmp_dir+head_file)
            except OSError:
                pass
    else:
        sys.stderr.write('Error: Please generate a fasta file with chromosome headers and either all chr prefixes or no chr prefixes'+'\n')
        if not debugging:
            shutil.rmtree(tmp_dir, ignore_errors = True)
            try:
                os.rmdir(out_dir)
            except OSError:
                pass
        sys.exit(1)
    return new_fasta_file

def get_sequence(local_fasta_file, instances_file, seq_instances_file, tmp_dir, new_fasta_file):
    if not new_fasta_file:
        if verbose:
            time_now('Extracting DNA sequence features..', starting = True)  
        os.system('(bedtools getfasta -fi %s -bed %s -name -bedOut -s > %s) 2>/dev/null' % (local_fasta_file, tmp_dir+instances_file, tmp_dir+seq_instances_file))
        if not os.path.isfile(tmp_dir+seq_instances_file):
            sys.stderr.write('Error: Unable to extract DNA sequences, perhaps something is wrong with the fasta file?'+'\n')
            if not debugging:
                shutil.rmtree(tmp_dir, ignore_errors = True)
                try:
                    os.rmdir(out_dir)
                except OSError:
                    pass
            sys.exit(1)
        if not debugging:
            try:
                os.remove(tmp_dir+instances_file) 
            except OSError:
                pass
        if not machine_learning_mode:
            try:
               os.remove(local_fasta_file+'.fai') 
            except OSError:
                pass
    else:
        if verbose:
            time_now('Extracting DNA sequence features..', starting = True)  
        os.system('(bedtools getfasta -fi %s -bed %s -name -bedOut -s > %s) 2>/dev/null' % (new_fasta_file, tmp_dir+instances_file, tmp_dir+seq_instances_file))
        if not os.path.isfile(tmp_dir+seq_instances_file):
            sys.stderr.write('Error: Unable to extract DNA sequences, perhaps something is wrong with the fasta file?'+'\n')
            if not debugging:
                shutil.rmtree(tmp_dir, ignore_errors = True)
                try:
                    os.rmdir(out_dir)
                except OSError:
                    pass
            sys.exit(1)
        if not debugging:
            try:
                os.remove(tmp_dir+instances_file) 
            except OSError:
                pass
        if not machine_learning_mode:
            try:
               os.remove(new_fasta_file+'.fai') 
            except OSError:
                pass
        if not machine_learning_mode and not debugging:
            try:
               os.remove(new_fasta_file) 
            except OSError:
                pass
        
def load_data_instances(window_merged_refined_utr_file_with_sequence, tmp_dir):
    utr = {}
    with open(tmp_dir+window_merged_refined_utr_file_with_sequence, 'r') as f:
        for line in f:
            if line.strip().split('\t')[0] == '\n' or line[0] == '#' or not line.strip():
                continue
            fields = line.strip('\n').split('\t')
            tpm = float(fields[4])
            bases = fields[8]
            split_name = fields[3].split('/')
            name = ('/').join(['%s|%s-%s' % (x, fields[6], fields[7]) for x in split_name])
            utr[name] = [fields[0], fields[5], fields[6], fields[7], tpm, bases]
    if not debugging:
        try:
            os.remove(tmp_dir+window_merged_refined_utr_file_with_sequence)
        except OSError:
            pass
    return utr
 
def extract_dna_features(utr):
    def signal_variant_indicator(bas, signal_variants):
        bas = bas.upper()
        in_string = [('').join(list(x)) for x in zip(bas, bas[1:], bas[2:], bas[3:], bas[4:], bas[5:])]
        if any(sig in signal_variants for sig in in_string):
            return 1
        else:
            return -1
    fl_dict = {}
    for ev in utr:
        if ev not in fl_dict:
            fl_dict[ev] = []
        curr = utr[ev]
        bases = curr[-1].upper()
        bases = bases[dna_probe_region_size-signal_window_size:dna_probe_region_size+window_size-min_upstream_signal_size]
        fl_dict[ev].append(signal_variant_indicator(bases, [canon_sigs[0]]))
        fl_dict[ev].append(signal_variant_indicator(bases, [canon_sigs[1]]))
        fl_dict[ev].append(signal_variant_indicator(bases, [canon_sigs[2]]))
        fl_dict[ev].append(signal_variant_indicator(bases, weak_sigs))                  
    return fl_dict
 
def base_freq_features(utr):                   
    def get_sig_loc(bas, signal_variants):
        in_string = [('').join(list(x)) for x in zip(bas, bas[1:], bas[2:], bas[3:], bas[4:], bas[5:])]
        for i in signal_variants:
            try:
                idx = in_string.index(i) 
                if idx:
                    return idx
            except ValueError:
                pass
        else:
            return None
    def get_string_gs(bas):
        in_string = [x for x in zip(bas, bas[1:], bas[2:], bas[3:], bas[4:], bas[5:])]
        count = 0
        for l in in_string:
            if l.count('G') >= 5:
                count += 1
        return count, len(in_string)
    def get_string_us(bas):
        in_string = [x for x in zip(bas, bas[1:], bas[2:])]
        count = 0
        for l in in_string:
            if l == ('T', 'T', 'T') or l == ('T', 'T', 'T'):
                count += 1
        return count, len(in_string)
    def get_gts(bas):
        in_string = [x for x in zip(bas, bas[1:])]
        count = 0
        for l in in_string:
            if l == ('G', 'T') or l == ('T', 'G'):
                count += 1
        return count, len(in_string)
    def get_gtgts(bas):
        in_string = [x for x in zip(bas, bas[1:], bas[2:], bas[3:])]
        count = 0
        for l in in_string:
            if l == ('G', 'T', 'G', 'T') or l == ('T', 'G', 'T', 'G'):
                count += 1
        return count, len(in_string)
    def get_aus(bas):
        in_string = [x for x in zip(bas, bas[1:])]
        count = 0
        for l in in_string:
            if l == ('A', 'T'):
                count += 1
        return count, len(in_string)  
    def get_us(bas):
        return bas.count('T'), len(bas)
    def get_tgta(bas):
        in_string = [x for x in zip(bas, bas[1:], bas[2:], bas[3:])]
        count = 0
        for l in in_string:
            if l == ('T', 'G', 'T', 'A') or l == ('T', 'A', 'T', 'A'):
                count += 1
        return count, len(in_string)
    base_combos = ['A', 'T', 'C', 'G']
    fl_dict = {}
    for ev in utr:
        if ev not in fl_dict:
            fl_dict[ev] = []
        curr = utr[ev]
        bases = curr[-1].upper()
        fin_bases = []
        for i in bases:
            if i not in base_combos:
                fin_bases.extend('N')
            else:
                fin_bases.extend(i)  
        sig_loc = get_sig_loc(fin_bases[dna_probe_region_size-signal_window_size:dna_probe_region_size+window_size-min_upstream_signal_size], canon_sigs)
        if sig_loc:
            sig_pos = sig_loc+dna_probe_region_size-signal_window_size
            up_us_count, up_us_len = get_us(fin_bases[sig_pos-up_us_adj:sig_pos])
            mid_us_count, mid_us_len = get_us(fin_bases[sig_pos+6:sig_pos+signal_window_size])
            gtgts_count, gtgts_len = get_gtgts(fin_bases[sig_pos+min_upstream_signal_size+6:sig_pos+signal_window_size+down_pos_adj])
            us_count, us_len = get_string_us(fin_bases[sig_pos+min_upstream_signal_size+6:sig_pos+signal_window_size+down_pos_adj])
            gts_count, gts_len = get_gts(fin_bases[sig_pos+min_upstream_signal_size+6:sig_pos+signal_window_size+down_pos_adj])
            tgta_count, tgta_len = get_tgta(fin_bases[sig_pos-tgtg_adj:sig_pos])
            if sig_pos+signal_window_size+100+6 <= len(fin_bases) and sig_pos+min_upstream_signal_size-100 >= 0:
                aus_count, aus_len = get_aus(fin_bases[sig_pos+min_upstream_signal_size-100:sig_pos+signal_window_size+100+6])
            elif sig_pos+signal_window_size+100+6 > len(fin_bases) and sig_pos+min_upstream_signal_size-100 >= 0:
                aus_count, aus_len = get_aus(fin_bases[sig_pos+min_upstream_signal_size-100:])
            elif sig_pos+signal_window_size+100+6 <= len(fin_bases) and sig_pos+min_upstream_signal_size-100 < 0:
                aus_count, aus_len = get_aus(fin_bases[:sig_pos+signal_window_size+100+6])
            else:
                aus_count, aus_len = get_aus(fin_bases)
            if sig_pos+min_upstream_signal_size+6+pos_gs_adj+window_size <= len(fin_bases):
                gs_count, gs_len = get_string_gs(fin_bases[sig_pos+min_upstream_signal_size+6+pos_gs_adj:])
            else:
                gs_count, gs_len = get_string_gs(fin_bases[sig_pos+min_upstream_signal_size+6+pos_gs_adj:sig_pos+min_upstream_signal_size+6+pos_gs_adj+window_size])  
        else:
            gtgts_count, gtgts_len = get_gtgts(fin_bases[dna_probe_region_size+up_pos_adj:dna_probe_region_size+window_size+down_pos_adj])
            us_count, us_len = get_string_us(fin_bases[dna_probe_region_size+up_pos_adj:dna_probe_region_size+window_size+down_pos_adj])
            gts_count, gts_len = get_gts(fin_bases[dna_probe_region_size+up_pos_adj:dna_probe_region_size+window_size+down_pos_adj])
            gs_count, gs_len = get_string_gs(fin_bases[dna_probe_region_size+pos_gs_adj:])
            aus_count, aus_len = get_aus(fin_bases)
            up_us_count, up_us_len = get_us(fin_bases[dna_probe_region_size-signal_window_size-up_us_adj:dna_probe_region_size+window_size-min_upstream_signal_size])
            mid_us_count, mid_us_len = get_us(fin_bases[dna_probe_region_size-signal_window_size:dna_probe_region_size+window_size])
            tgta_count, tgta_len = get_tgta(fin_bases[dna_probe_region_size-signal_window_size-tgtg_adj:dna_probe_region_size+window_size-min_upstream_signal_size])
        if gtgts_count/gtgts_len > gtgts_cut:
            fl_dict[ev].append(1)
        else:
            fl_dict[ev].append(-1)
        if gts_count/gts_len > gts_cut:
            fl_dict[ev].append(1)
        else:
            fl_dict[ev].append(-1)
        if us_count/us_len > us_cut:
            fl_dict[ev].append(1)
        else:
            fl_dict[ev].append(-1)   
        if gs_count/gs_len > gs_cut:
            fl_dict[ev].append(1)
        else:
            fl_dict[ev].append(-1)  
        if aus_count/aus_len > aus_cut:
            fl_dict[ev].append(1)
        else:
            fl_dict[ev].append(-1)  
        if up_us_count/up_us_len > up_us_cut:
            fl_dict[ev].append(1)
        else:
            fl_dict[ev].append(-1)  
        if mid_us_count/mid_us_len > mid_us_cut:
            fl_dict[ev].append(1)
        else:
            fl_dict[ev].append(-1)  
        if tgta_count/tgta_len > tgta_cut:
            fl_dict[ev].append(1)
        else:
            fl_dict[ev].append(-1) 
    if verbose:
        time_now('', finishing = True)
    return fl_dict
 
def extract_rna_seq_coverage(bedgraph_files, tmp_dir, utr):
    if verbose:
        time_now('Extracting RNA sequencing features..', starting = True)  
    fl_dict = {}
    All_Samples_Total_depth = []
    cur_sample_total_depth = 0
    curr_sample_All_chroms_coverage_dict = {}
    curr_sample_All_chroms_coverage_dict_minus = {}
    with open(tmp_dir+bedgraph_files, 'r') as f:
        for line in f:
            if line.strip().split('\t')[0] == '\n' or line[0] == '#' or not line.strip():
                continue
            fields = line.strip('\n').split('\t')
            chrom_name = fields[0]
            if not chrom_name.startswith('chr'):
                chrom_name = 'chr' + chrom_name
            region_start = int(float(fields[1]))
            region_end = int(float(fields[2]))
            cur_sample_total_depth += int(float(fields[-2])) * (region_end - region_start) + int(float(fields[-1])) * (region_end - region_start)
            if chrom_name not in curr_sample_All_chroms_coverage_dict:
                curr_sample_All_chroms_coverage_dict[chrom_name] = [[0],[0]]
            if region_start > curr_sample_All_chroms_coverage_dict[chrom_name][0][-1]:
                curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_start)
                curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)
            curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_end)
            curr_sample_All_chroms_coverage_dict[chrom_name][1].append(int(float(fields[-2])))
            if chrom_name not in curr_sample_All_chroms_coverage_dict_minus:
                curr_sample_All_chroms_coverage_dict_minus[chrom_name] = [[0],[0]]
            if region_start > curr_sample_All_chroms_coverage_dict_minus[chrom_name][0][-1]:
                curr_sample_All_chroms_coverage_dict_minus[chrom_name][0].append(region_start)
                curr_sample_All_chroms_coverage_dict_minus[chrom_name][1].append(0)
            curr_sample_All_chroms_coverage_dict_minus[chrom_name][0].append(region_end)
            curr_sample_All_chroms_coverage_dict_minus[chrom_name][1].append(int(float(fields[-1])))
        curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)
        curr_sample_All_chroms_coverage_dict_minus[chrom_name][1].append(0)
    if cur_sample_total_depth != 0:
        All_Samples_Total_depth.append(cur_sample_total_depth)
    else:
        sys.stderr.write('Error: No read coverage detected'+'\n')
        if not debugging:
            shutil.rmtree(tmp_dir, ignore_errors = True)
            try:
                os.rmdir(out_dir)
            except OSError:
                pass
        sys.exit(1)
    for ev in utr:
        if ev not in fl_dict:
            fl_dict[ev] = []
        curr = utr[ev]
        curr_chr = curr[0]
        strand = curr[1]
        inst_start = int(float(curr[2]))
        inst_end = int(float(curr[3]))
        region_start = int(float(inst_start - rna_seq_probe_region_size))
        region_end = int(float(inst_end + rna_seq_probe_region_size))
        if strand == '+' and curr_chr in curr_sample_All_chroms_coverage_dict:
            curr_chr_coverage = curr_sample_All_chroms_coverage_dict[curr_chr]
            left_region_index = bisect(curr_chr_coverage[0], region_start)
            right_region_index = bisect(curr_chr_coverage[0], region_end)
            extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]
            extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
            extracted_3UTR_region.insert(0, region_start)
            extracted_3UTR_region.append(region_end)
            bp_coverage = np.zeros(extracted_3UTR_region[-1] - extracted_3UTR_region[0])
        elif strand == '-' and curr_chr in curr_sample_All_chroms_coverage_dict_minus:
            curr_chr_coverage = curr_sample_All_chroms_coverage_dict_minus[curr_chr]
            left_region_index = bisect(curr_chr_coverage[0], region_start)
            right_region_index = bisect(curr_chr_coverage[0], region_end)
            extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]
            extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
            extracted_3UTR_region.insert(0, region_start)
            extracted_3UTR_region.append(region_end)
            bp_coverage = np.zeros(extracted_3UTR_region[-1] - extracted_3UTR_region[0])
        else:
            continue
        relative_start = extracted_3UTR_region[0]
        for i in range(len(extracted_coverage)):
            curr_region_start = extracted_3UTR_region[i] - relative_start
            curr_region_end = extracted_3UTR_region[i+1] - relative_start
            bp_coverage[curr_region_start:curr_region_end] = extracted_coverage[i]
        if strand == '-':
            bp_coverage = bp_coverage[::-1]
        fl_dict[ev].append(bp_coverage)
    if not machine_learning_mode and not debugging:
        try:
            os.remove(tmp_dir+bedgraph_files) 
        except OSError:
            pass
    return fl_dict, All_Samples_Total_depth
 
def extract_local_rna_seq_features(rna_seq_cov_dict, All_Samples_Total_depth):
    fl_dict = {}
    split_sizes = [window_size/num_split]*num_split
    split_sizes[0] = np.ceil(split_sizes[0])
    split_sizes = [int(float(np.floor(x))) for x in split_sizes]
    for ev in rna_seq_cov_dict:
        exp_list = []
        cov_array = rna_seq_cov_dict[ev][0]
        reg_one = cov_array[rna_seq_probe_region_size:rna_seq_probe_region_size+split_sizes[0]]
        reg_two = cov_array[rna_seq_probe_region_size+split_sizes[0]:rna_seq_probe_region_size+split_sizes[0]+split_sizes[1]]
        reg_three = cov_array[rna_seq_probe_region_size+split_sizes[0]+split_sizes[1]:rna_seq_probe_region_size+split_sizes[0]+split_sizes[1]+split_sizes[2]]
        reg_one_med = np.median(reg_one)
        reg_two_med = np.median(reg_two)
        reg_three_med = np.median(reg_three)
        diff_one_two = reg_one_med - reg_two_med
        diff_two_three = reg_two_med - reg_three_med
        if reg_one_med+reg_two_med+reg_three_med != 0:
            rat_reg_one = reg_one_med/(reg_one_med+reg_two_med+reg_three_med)
            rat_reg_two = reg_two_med/(reg_one_med+reg_two_med+reg_three_med)
            rat_reg_three = reg_three_med/(reg_one_med+reg_two_med+reg_three_med)
        else:
            rat_reg_one = 0
            rat_reg_two = 0
            rat_reg_three = 0
        if reg_one_med+reg_three_med != 0:
            no_mid_rat_reg_one = reg_one_med/(reg_one_med+reg_three_med)
            no_mid_rat_reg_three = reg_three_med/(reg_one_med+reg_three_med)
        else:
            no_mid_rat_reg_one = 0
            no_mid_rat_reg_three = 0
        exp_list.extend([diff_one_two] + [diff_two_three] + [rat_reg_one] + [rat_reg_two] + [rat_reg_three] + [no_mid_rat_reg_one] + [no_mid_rat_reg_three])
        fl_dict[ev] = exp_list
    return fl_dict
  
def extract_rna_seq_features(rna_seq_cov_dict, All_Samples_Total_depth):
    fl_dict = {}
    for ev in rna_seq_cov_dict:
        exp_list = []
        cov_array = rna_seq_cov_dict[ev][0]
        reg_one_med = np.median(cov_array[0:rna_seq_probe_region_size])
        reg_two_med = np.median(cov_array[rna_seq_probe_region_size:rna_seq_probe_region_size+window_size])
        reg_three_med = np.median(cov_array[rna_seq_probe_region_size+window_size:rna_seq_probe_region_size*2+window_size])
        diff_one_two = reg_one_med - reg_two_med
        diff_two_three = reg_two_med - reg_three_med
        if reg_one_med+reg_two_med+reg_three_med != 0:
            rat_reg_one = reg_one_med/(reg_one_med+reg_two_med+reg_three_med)
            rat_reg_two = reg_two_med/(reg_one_med+reg_two_med+reg_three_med)
            rat_reg_three = reg_three_med/(reg_one_med+reg_two_med+reg_three_med)
        else:
            rat_reg_one = 0
            rat_reg_two = 0
            rat_reg_three = 0
        if reg_one_med+reg_three_med != 0:
            no_mid_rat_reg_one = reg_one_med/(reg_one_med+reg_three_med)
            no_mid_rat_reg_three = reg_three_med/(reg_one_med+reg_three_med)
        else:
            no_mid_rat_reg_one = 0
            no_mid_rat_reg_three = 0
        exp_list.extend([diff_one_two] + [diff_two_three] + [rat_reg_one] + [rat_reg_two] + [rat_reg_three] + [no_mid_rat_reg_one] + [no_mid_rat_reg_three])
        fl_dict[ev] = exp_list
    if verbose:
        time_now('', finishing = True)
    return fl_dict

def stringtie_ends_as_dict(b_list):
    if verbose:
        time_now('Extracting reconstruction end locations feature..', starting = True)  
    label_dict = {}
    label_dict_minus = {}
    label_dict_fin = {}
    label_dict_minus_fin = {}
    for i in b_list:
        chrom_name = i[0]
        strand = i[5]
        if strand == '+':
            polya_site = int(float(i[2]))
        elif strand == '-':
            polya_site = int(float(i[1]))
        name = str(polya_site)+chrom_name
        if not chrom_name.startswith('chr'):
            chrom_name = 'chr' + chrom_name
        if strand == '-':
            if name not in label_dict_minus:
                label_dict_minus[name] = []
                label_dict_minus[name].extend([polya_site, chrom_name])
        elif strand == '+':
            if name not in label_dict:
                label_dict[name] = []
                label_dict[name].extend([polya_site, chrom_name])
    for i in label_dict:
        if label_dict[i][1] not in label_dict_fin:
            label_dict_fin[label_dict[i][1]] = []
        label_dict_fin[label_dict[i][1]].extend([label_dict[i][0]])
    for i in label_dict_minus:
        if label_dict_minus[i][1] not in label_dict_minus_fin:
            label_dict_minus_fin[label_dict_minus[i][1]] = []
        label_dict_minus_fin[label_dict_minus[i][1]].extend([label_dict_minus[i][0]])
    return label_dict_fin, label_dict_minus_fin
 
def extract_stringtie_ends_as_feature(gencode_label_dict, gencode_label_dict_minus, utr):
    fl_dict = {}
    for ev in utr:
        if ev not in fl_dict:
            fl_dict[ev] = []
        curr = utr[ev]
        curr_chr = curr[0]
        strand = curr[1]
        region_start = curr[2]
        region_end = curr[3]
        if strand == '+':
            if curr_chr in gencode_label_dict:
                probe_start, probe_end = int(float(region_start)), int(float(region_end))
                if any([i for i in gencode_label_dict[curr_chr] if probe_start <= i < probe_end]):
                    fl_dict[ev].append(1)
                else:
                    fl_dict[ev].append(-1)
            else:
                fl_dict[ev].append(-1)
        elif strand == '-':
            if curr_chr in gencode_label_dict_minus:
                probe_start, probe_end = int(float(region_start)), int(float(region_end))
                if any([i for i in gencode_label_dict_minus[curr_chr] if probe_start <= i < probe_end]):
                    fl_dict[ev].append(1)
                else:
                    fl_dict[ev].append(-1)
            else:
                fl_dict[ev].append(-1)
    if verbose:
        time_now('', finishing = True)
    return fl_dict

def extract_label_locations(annotation_file, tmp_dir):
    if machine_learning_mode:
        if verbose:
            time_now('Machine learning mode enabled, loading polyA sites and labeling data instances..', starting = True)    
        label_dict = {}
        label_dict_minus = {}
        try:
            with open(annotation_file, 'r') as f:
                for line in f:
                    if line.strip().split('\t')[0] == '\n' or line[0] == '#' or not line.strip():
                        continue
                    fields = line.strip().split('\t')
                    chrom_name = fields[chrom_loc]
                    strand = fields[strand_loc]
                    polya_site = int(float(fields[site_loc]))
                    if not chrom_name.startswith('chr'):
                        chrom_name = 'chr' + chrom_name
                    if strand == '-':
                        if chrom_name not in label_dict_minus:
                            label_dict_minus[chrom_name] = [polya_site]
                        else:
                            label_dict_minus[chrom_name].extend([polya_site])
                    elif strand == '+':
                        if chrom_name not in label_dict:
                            label_dict[chrom_name] = [polya_site]
                        else:
                            label_dict[chrom_name].extend([polya_site])
            for key in label_dict:
                label_dict[key] = list(set(label_dict[key]))
                label_dict[key].sort()
            for key in label_dict_minus:
                label_dict_minus[key] = list(set(label_dict_minus[key]))
                label_dict_minus[key].sort()
            return label_dict, label_dict_minus
        except:
            sys.stderr.write('Error: Could not open label file or could not properly read in data'+'\n')
            if not debugging:
                shutil.rmtree(tmp_dir, ignore_errors = True)
                try:
                    os.rmdir(out_dir)
                except OSError:
                    pass
            sys.exit(1)
        if len(label_dict) + len(label_dict_minus) == 0:
            sys.stderr.write('Error: label info is empty: Perhaps the file is not tab separated and/or the columns denoting polyA sites, chromosome, and strand are wrong?'+'\n')
            if not debugging:
                shutil.rmtree(tmp_dir, ignore_errors = True)
                try:
                    os.rmdir(out_dir)
                except OSError:
                    pass
            sys.exit(1)
    else:
        return None, None

def extract_labels(gencode_label_dict, gencode_label_dict_minus, utr):
    if machine_learning_mode:
        fl_dict = {}
        for ev in utr:
            if ev not in fl_dict:
                fl_dict[ev] = []
            curr = utr[ev]
            curr_chr = curr[0]
            strand = curr[1]
            region_start = curr[2]
            region_end = curr[3]
            if strand == '+':
                if curr_chr in gencode_label_dict:
                    probe_start, probe_end = int(float(region_start)), int(float(region_end))
                    if any([i for i in gencode_label_dict[curr_chr] if probe_start <= i < probe_end]):
                        fl_dict[ev].append(1)
                    else:
                        fl_dict[ev].append(0)
                else:
                    fl_dict[ev].append(0)
            elif strand == '-':
                if curr_chr in gencode_label_dict_minus:
                    probe_start, probe_end = int(float(region_start)), int(float(region_end))
                    if any([i for i in gencode_label_dict_minus[curr_chr] if probe_start <= i < probe_end]):
                        fl_dict[ev].append(1)
                    else:
                        fl_dict[ev].append(0)
                else:
                    fl_dict[ev].append(0)
        if verbose:
            time_now('', finishing = True)
        return fl_dict
    else:
        return None

def combine_feat_lab_dicts(*args):
    result = {}
    for dic in args:
        if not dic == None:
            for key in (result.keys() | dic.keys()):
                if key in dic:
                    if type(dic[key]) is list:
                        result.setdefault(key, []).extend(dic[key])
                    else:
                        result.setdefault(key, []).append(dic[key])
    return result

def bi_di_format(in_dict, col_names):
    out_frame = pd.DataFrame.from_dict(in_dict, orient = 'index')
    out_frame.columns = col_names
    trans_name = []
    strand = []
    pos = []
    for i in out_frame.index:
        idx_name = i
        long_name = idx_name.split('/')
        trans_name.append([y.rsplit('|', 1)[0] for y in long_name])
        strand.append(idx_name.split('/')[0].rsplit('|', 2)[1])
        pos.append(idx_name.split('|')[-1])
    out_frame['trans'] = [(',').join(x) for x in trans_name]
    out_frame['strand'] = strand
    out_frame['pos'] = pos
    out_frame['index'] = out_frame.index
    pos_out_frame = out_frame[out_frame['strand'] == '+']
    neg_out_frame = out_frame[out_frame['strand'] == '-']
    pos_out_frame = pos_out_frame.sort_values('index', ascending = True)
    neg_out_frame = neg_out_frame.sort_values('index', ascending = False)
    out_frame = pd.concat([pos_out_frame, neg_out_frame])
    out_frame.drop('index', axis = 1, inplace = True)
    out_frame.drop('strand', axis = 1, inplace = True)
    out_frame = out_frame.groupby('trans', as_index = False).agg(lambda x: list(x))
    out_list = out_frame.values.tolist()
    return out_list

def split_names_data(in_list):
    fin_list = []
    fin_names = []
    for i in range(len(in_list)):
        fin_list.append(np.transpose(np.array(in_list[i][1:-1])))
        fin_names.append([in_list[i][0], in_list[i][-1]])
    return fin_list, fin_names

def get_continuous_vars(in_list):
    all_idx = np.arange(np.array(in_list)[-1].shape[-1])
    con_idx = []
    for i in range(12, 26):
        con_idx.append(i)
    bi_idx = []
    for x in all_idx:
        if x not in con_idx:
            bi_idx.append(x)
    return con_idx, bi_idx

def bi_di_format_features(in_list, con_idx, pad_val):
    def padding(train_X, test_X, train_y, test_y, val, max_length):  
        padded_train_X = pad_sequences(train_X, padding = 'post', value = val, dtype = 'float64', maxlen = max_length)
        padded_test_X = pad_sequences(test_X, padding = 'post', value = val, dtype = 'float64', maxlen = max_length)
        padded_train_y = pad_sequences(train_y, padding = 'post', value = val, dtype = 'float64', maxlen = max_length)
        padded_test_y = pad_sequences(test_y, padding = 'post', value = val, dtype = 'float64', maxlen = max_length)
        return padded_train_X, padded_test_X, padded_train_y, padded_test_y
    x_vals = []
    y_vals = []
    for i in range(len(in_list)): 
        x_vals.append([x[:-1] for x in in_list[i]])
        y_vals.append([x[-1] for x in in_list[i]])
    fin_x_vals = [] 
    fin_y_vals = []
    for a, b in zip(range(len(x_vals)), range(len(y_vals))):
        fin_x_vals.append(np.concatenate([x_vals[a]]))
        fin_y_vals.append(np.concatenate([y_vals[b]]))
    if not rand_seed:
        X_tr_bi, X_te_bi, y_tr_bi, y_te_bi, idx_tr, idx_te = train_test_split(fin_x_vals, fin_y_vals, range(len(fin_x_vals)), test_size = 0.2)
    else:
        X_tr_bi, X_te_bi, y_tr_bi, y_te_bi, idx_tr, idx_te = train_test_split(fin_x_vals, fin_y_vals, range(len(fin_x_vals)), random_state = rand_seed, test_size = 0.2) 
    scalers = {}   
    for c in con_idx: 
        scalers[c] = StandardScaler() 
        tmp = []     
        for d in range(len(X_tr_bi)):
            for e in range(len(X_tr_bi[d])):
                tmp.append(X_tr_bi[d][e][c])
        scalers[c] = scalers[c].fit(np.array(tmp).reshape(-1, 1))
        for f in range(len(X_tr_bi)):
            for g in range(len(X_tr_bi[f])):
                X_tr_bi[f][g][c] = scalers[c].transform(X_tr_bi[f][g][c].reshape(-1, 1))
        for h in range(len(X_te_bi)):
            for j in range(len(X_te_bi[h])):
                X_te_bi[h][j][c] = scalers[c].transform(X_te_bi[h][j][c].reshape(-1, 1))            
    if tmp_maxlen == 'default':
        max_len = []
        for k in range(len(X_tr_bi)):
            max_len.append(len(X_tr_bi[k])) 
        for l in range(len(X_te_bi)):  
            max_len.append(len(X_te_bi[l]))    
        fin_maxlen = max(max_len) 
    else:
        fin_maxlen = tmp_maxlen
    fin_X_tr_bi, fin_X_te_bi, fin_y_tr_bi, fin_y_te_bi = padding(X_tr_bi, X_te_bi, y_tr_bi, y_te_bi, pad_val, fin_maxlen)
    fin_y_tr_bi = np.expand_dims(fin_y_tr_bi, axis = 2)
    fin_y_te_bi = np.expand_dims(fin_y_te_bi, axis = 2)
    return fin_X_tr_bi, fin_X_te_bi, fin_y_tr_bi, fin_y_te_bi, idx_tr, idx_te, fin_X_tr_bi.shape[-1], fin_maxlen, scalers

def generate_weights(y_tr_bi, pad_val):
    return dict(enumerate(class_weight.compute_class_weight(class_weight = 'balanced', classes = np.unique(y_tr_bi)[np.unique(y_tr_bi) != pad_val], y = y_tr_bi.flatten()[y_tr_bi.flatten() != pad_val])))

def dist_class_weights(weights_dict, y_tr_bi):
    dist_y_tr_bi = y_tr_bi.copy()
    fin_dist_y_tr_bi = np.array(list(map(np.vectorize(weights_dict.get), dist_y_tr_bi)))
    return fin_dist_y_tr_bi

def save_files(file_list, name_list, tmp_dir):
    for a, b in zip(file_list, name_list):
        file = open(tmp_dir+b, 'wb')
        pickle.dump(a, file)
        file.close()
        
def open_files(name_list):
    file = open(name_list, 'rb')
    fi = pickle.load(file)
    file.close()
    return fi

def format_features(in_list, con_idx, pad_val, scale_dict, max_len):
    def padding(train_X, val, max_length):  
        padded_train_X = pad_sequences(train_X, padding = 'post', value = val, dtype = 'float64', maxlen = max_length)
        return padded_train_X
    X_tr_bi = [] 
    if machine_learning_mode:
        x_vals = []
        for i in range(len(in_list)): 
            x_vals.append([x[:-1] for x in in_list[i]])
        for a in range(len(x_vals)):
            X_tr_bi.append(np.concatenate([x_vals[a]]))
    else:
        for a in range(len(in_list)):
            X_tr_bi.append(np.concatenate([in_list[a]]))
    for c in con_idx: 
        for f in range(len(X_tr_bi)):
            for g in range(len(X_tr_bi[f])):
                X_tr_bi[f][g][c] = scale_dict[c].transform(X_tr_bi[f][g][c].reshape(-1, 1))  
    fin_X_tr_bi = padding(X_tr_bi, pad_val, max_len)
    return fin_X_tr_bi

def keras_load_mod(in_name):
    return load_model(in_name)

def new_define_data_instances(tmp_dir, instances_file, chrom_sizes_dict):
    if verbose:
        time_now('Redefining data instances..', starting = True)  
    window_refined_utrs_file = 'window.bed'
    output_write = open(tmp_dir+instances_file, 'w')
    with open(tmp_dir+window_refined_utrs_file, 'r') as f:
        for line in f:
            if line.strip().split('\t')[0] == '\n' or line[0] == '#' or not line.strip():
                continue
            fields = line.strip('\n').split('\t')
            name = fields[3]
            start = int(float(fields[1]))
            end = int(float(fields[2]))
            strand = fields[5]
            tpm = fields[4]
            comb_tpm = str(sum([float(x) for x in tpm.split('/')]))
            for i in range(start, end, window_size):
                base_start = i - dna_probe_region_size
                base_end = i + dna_probe_region_size + window_size
                if base_start < 0:
                    continue
                if fields[0] in chrom_sizes_dict:
                    total_chrom_length = chrom_sizes_dict[fields[0]][0]
                    if int(float(total_chrom_length)) < base_end:
                        continue
                write_line = [fields[0], str(base_start), str(base_end), name, comb_tpm, strand, str(i), str(i + window_size)]
                output_write.writelines('\t'.join(write_line) + '\n')
    output_write.close()
    if not debugging:
        try:
            os.remove(tmp_dir+window_refined_utrs_file)
        except OSError:
            pass
    if verbose:
        time_now('', finishing = True) 
 
def trunc_names(in_names, length):
    fin_list = []
    for i in range(len(in_names)):
        tmp_names = in_names[i][1][0:length]
        fin_list.append([in_names[i][0], tmp_names])
    return fin_list

def add_prob_pred(test_X, length, in_mod):
    pred = in_mod.predict_classes(test_X)
    yhat = in_mod.predict_proba(test_X)
    fin_X = np.zeros((len(test_X), length, test_X.shape[-1]+2))
    for a, b, c, d in zip(yhat, pred, test_X, range(len(fin_X))):
        tmp = np.concatenate([c, a, b], axis = 1)
        fin_X[d] = tmp
    return fin_X

def remove_mask_add_names(test_X, pad_val, in_names):
    fin_list = []
    idxs = []
    for i in range(len(test_X)):
        idxs.append(np.array([x for x in range(len(np.where(test_X[i, :, 0] != pad_val)[0]))]))
    for a, b in zip(range(len(test_X)), idxs):
        fin_list.append([np.concatenate([test_X[a, b, -2:], np.array(in_names[a][1]).reshape(-1, 1)], axis = 1), in_names[a][0]])
    return fin_list

def final_gtf(in_preds, gtf_file, fo, thresh, tmp_dir):
    def bed_line(estart_val, eend_val, field_vals, nline_val, in_att, in_lab):
        estp = estart_val[0] - 1
        eedp = eend_val[-1]
        curr_strand = field_vals[6]
        curr_chrom = field_vals[0]
        out_line = [curr_chrom, str(estp), str(eedp), str(nline_val), str(in_att), curr_strand, str(estp), str(eedp), str(in_lab), str(len(estart_val))] 
        seglen = [eend_val[i] - estart_val[i] + 1 for i in range(len(estart_val))]
        segstart = [estart_val[i] - estart_val[0] for i in range(len(estart_val))]
        strl = str(seglen[0])
        for i in range(1, len(seglen)):
            strl += ',' + str(seglen[i])
        strs = str(segstart[0])
        for i in range(1, len(segstart)):
            strs +=',' + str(segstart[i])
        out_line.extend([strl, strs])
        return out_line
    def write_gtf_lines(elems, in_id, num_trans):
        fin_list = []                   
        chrom = elems[0]
        att = elems[4]
        source_name = elems[8]
        attribs = att.split(',')
        sources = source_name.split(',')
        start = int(float(elems[1])) + 1
        end = int(float(elems[2]))
        strand = elems[5]
        ver = int(float(elems[3]))
        fin_attribs = []
        if 'aptardi' in sources:
            for l in attribs:
                new_name = re.findall(r'transcript_id \"([\w\.]+)\"', l)
                if len(new_name) != 0:
                    new_name.extend(['.', str(ver)])
                    new_name = ('').join(new_name)
                    fin_attribs.append(re.sub(r'transcript_id \"([\w\.]+)\"', 'transcript_id '+'"{}"'.format(new_name), l))
                else:
                    new_name = ''.join(['aptardi.'+str(num_trans)+'.'+str(ver)])
                    fin_attribs.append(l+' transcript_id '+ '"{}"'.format(new_name))      
        else:
            for l in attribs:
                fin_attribs.append(l)
        block_count = int(elems[9])
        block_sizes = elems[10].split(',')
        block_starts = elems[11].split(',')
        count = 0
        fin_list.append([chrom, sources[0], 'transcript', str(start), str(end), '1000', strand, '.', fin_attribs[0]])
        for j in range(block_count):
            exon_start = int(float(start)) + int(float(block_starts[j]))
            exon_end = exon_start + int(float(block_sizes[j])) - 1
            group = fin_attribs[1:][j]
            source_group = sources[1:][j]
            count += 1
            fin_list.append([chrom, source_group, 'exon', str(exon_start), str(exon_end), '1000', strand, '.', group])
        return fin_list
    tot_line = 0
    estart = []
    eend = []
    prevfield = []
    all_att = []
    all_label = []
    prevtransid = ''
    prevstrand = ''
    prevgeneid = ''
    fo.writelines('# aptardi version %s\n' % (soft_ver))
    for lines in gtf_file:
        if lines[0] == '#':
            fo.writelines(lines)
            continue
        if lines.strip().split('\t')[0] == '\n' or not lines.strip():
            continue
        fields = lines.strip().split('\t')
        fields[8:] = [(' ').join(fields[8:])]
        att_field = fields[8]
        recon_field = fields[1]
        type_field = fields[2]
        strand = fields[6]
        transid = re.findall(r'transcript_id \"([\w\.]+)\"', att_field) 
        if len(transid) > 0:
            transid = transid[0]
        else:
            transid = ''
        geneid = re.findall(r'gene_id \"([\w\.]+)\"', att_field) 
        if len(geneid) > 0:
            geneid = geneid[0]
        else:
            geneid = ''
        if type_field == 'transcript' or (prevtransid != '' and transid != '' and transid != prevtransid):
            if len(estart) != 0:
                if len([i for i in in_preds if i[1].split('|')[0] == prevtransid and i[1].split('|')[1] == prevgeneid]) == 0:
                    tot_line += 1
                    line = bed_line(estart, eend, prevfield, 0, ','.join(all_att), ','.join(all_label))
                    all_lines = write_gtf_lines(line, prevtransid, tot_line)
                    for fin_lines in all_lines:
                        fo.writelines('\t'.join(fin_lines) + '\n')
                else:
                    curr_preds = [i for i in in_preds if i[1].split('|')[0] == prevtransid and i[1].split('|')[1] == prevgeneid][0]
                    tr_line = 0
                    tot_line += 1
                    line = bed_line(estart, eend, prevfield, 0, ','.join(all_att), ','.join(all_label))
                    all_lines = write_gtf_lines(line, prevtransid, tot_line)
                    for fin_lines in all_lines:
                        fo.writelines('\t'.join(fin_lines) + '\n')
                    for x in curr_preds[0]:
                        polya_prob = float(x[0])
                        if thresh == 'default':
                            polya_pred =  int(float(x[1]))
                        else:
                            if polya_prob >= thresh:
                                polya_pred = 1
                            else:
                                polya_pred = 0
                        if polya_pred != 1:
                            continue
                        reg_start, reg_end = [int(float(y)) for y in x[-1].split('-')]
                        if prevstrand == '-' and not reg_start+1 <= int(float(estart[0])) <= reg_end and not abs(reg_start+1 - int(float(estart[0]))) <= 100:
                            tr_line += 1
                            line = bed_line([reg_start+1]+estart[1:], eend, prevfield, tr_line, ','.join(all_att), ','.join(['aptardi']*2+all_label[2:]))
                            all_lines = write_gtf_lines(line, prevtransid, tot_line)
                            for fin_lines in all_lines:
                                fo.writelines('\t'.join(fin_lines) + '\n')
                        elif prevstrand == '+' and not reg_start+1 <= int(float(eend[-1])) <= reg_end and not abs(reg_end - int(float(eend[-1]))) <= 100:
                            tr_line += 1
                            line = bed_line(estart, eend[:-1]+[reg_end], prevfield, tr_line, ','.join(all_att), ','.join(['aptardi']+all_label[1:-1]+['aptardi']))
                            all_lines = write_gtf_lines(line, prevtransid, tot_line)
                            for fin_lines in all_lines:
                                fo.writelines('\t'.join(fin_lines) + '\n')
            estart = []
            eend = []
            all_att = []
            all_label = []
        prevfield = fields
        prevtransid = transid
        prevgeneid = geneid
        prevstrand = strand
        all_att.append(att_field)
        all_label.append(recon_field)
        if type_field == 'exon':
            est = int(float(fields[3]))
            eed = int(float(fields[4]))
            estart += [est]
            eend += [eed]
    if len(estart) != 0:
        if len([i for i in in_preds if i[1].split('|')[0] == prevtransid and i[1].split('|')[1] == prevgeneid]) == 0:
            tot_line += 1
            line = bed_line(estart, eend, prevfield, 0, ','.join(all_att), ','.join(all_label))
            all_lines = write_gtf_lines(line, prevtransid, tot_line)
            for fin_lines in all_lines:
                fo.writelines('\t'.join(fin_lines) + '\n')
        else:
            curr_preds = [i for i in in_preds if i[1].split('|')[0] == prevtransid and i[1].split('|')[1] == prevgeneid][0]
            tr_line = 0
            tot_line += 1
            line = bed_line(estart, eend, prevfield, 0, ','.join(all_att), ','.join(all_label))
            all_lines = write_gtf_lines(line, prevtransid, tot_line)
            for fin_lines in all_lines:
                fo.writelines('\t'.join(fin_lines) + '\n')
            for x in curr_preds[0]:
                polya_prob = float(x[0])
                if thresh == 'default':
                    polya_pred =  int(float(x[1]))
                else:
                    if polya_prob >= thresh:
                        polya_pred = 1
                    else:
                        polya_pred = 0
                if polya_pred != 1:
                    continue
                reg_start, reg_end = [int(float(y)) for y in x[-1].split('-')]
                if prevstrand == '-' and not reg_start+1 <= int(float(estart[0])) <= reg_end and not abs(reg_start+1 - int(float(estart[0]))) <= 100:
                    tr_line += 1
                    line = bed_line([reg_start+1]+estart[1:], eend, prevfield, tr_line, ','.join(all_att), ','.join(['aptardi']*2+all_label[2:]))
                    all_lines = write_gtf_lines(line, prevtransid, tot_line)
                    for fin_lines in all_lines:
                        fo.writelines('\t'.join(fin_lines) + '\n')
                elif prevstrand == '+' and not reg_start+1 <= int(float(eend[-1])) <= reg_end and not abs(reg_end - int(float(eend[-1]))) <= 100:
                    tr_line += 1
                    line = bed_line(estart, eend[:-1]+[reg_end], prevfield, tr_line, ','.join(all_att), ','.join(['aptardi']+all_label[1:-1]+['aptardi']))
                    all_lines = write_gtf_lines(line, prevtransid, tot_line)
                    for fin_lines in all_lines:
                        fo.writelines('\t'.join(fin_lines) + '\n')
    fo.close()
    if not debugging:
        shutil.rmtree(tmp_dir, ignore_errors = True)
        
utr_extension = 10000
start_utr_cutoff = 0.10
per_of_beg_utr = 0.05
per_bases_curr_window = 0.8
min_num_regions = 3 
rna_seq_probe_region_size = 100
num_split = 3
dna_probe_region_size = 100
signal_window_size = 36
min_upstream_signal_size = 7  
up_pos_adj = 10
down_pos_adj = 40
pos_gs_adj = 30
tgtg_adj = 40
up_us_adj = 50
gts_cut = 1/16*2*2
gtgts_cut = 1/256*2*6
us_cut = 1/64*8
aus_cut = 1/16*2
gs_cut = 24/4096*10
up_us_cut = 1/4*1.5
mid_us_cut = 1/4*1.5
tgta_cut = 1/256*2*6
canon_sigs = ['AATAAA', 'ATTAAA', 'AGTAAA']
weak_sigs = ['AAGAAA', 'AAAAAG', 'AATACA', 'TATAAA', 'ACTAAA', 'GATAAA', 'AATATA', 'CATAAA', 'AATAGA']
soft_ver = '1.4'
chr_sizes_name = 'chrom.sizes'
bg_file = 'rna_seq_cov.bg'     

def main():
    global orientation
    global debugging
    global machine_learning_mode
    global verbose
    global chrom_loc
    global strand_loc
    global site_loc
    global rand_seed
    global tmp_maxlen
    global window_size
    global out_dir
    args = parse_args() 
    orientation = args.a
    stg_gtf_wrap = args.r
    stg_gtf = stg_gtf_wrap.readlines()
    fasta_file = args.f
    bam_file = args.b
    debugging = args.d
    machine_learning_mode = args.m
    res_file = args.g
    threshold = args.p
    window_size = args.w
    rand_seed = args.c
    out_dir, temp_dir = args.o
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)
    else:
        print('Temporary directory name using date and time to seconds already exists, exiting')
        sys.exit(1)
    try:
        with open(temp_dir+'testing.txt', 'w') as f:
            f.write('testing')
    except OSError:
        print('Temporary directory is not writable, exiting')
        shutil.rmtree(temp_dir, ignore_errors = True)
        try:
            os.rmdir(out_dir)
        except OSError:
            pass
        sys.exit(1)
    try:
        with open(temp_dir+'testing.txt', 'r') as f:
            for line in f:
                test = line.strip().split('\t')
    except OSError:
        print('Temporary directory is not readable, exiting')
        shutil.rmtree(temp_dir, ignore_errors = True)
        try:
            os.rmdir(out_dir)
        except OSError:
            pass
        sys.exit(1)
    os.remove(temp_dir+'testing.txt')
    chrom_loc, strand_loc, site_loc, polya_sites, mod_name, sc_name, tmp_maxlen = check_ml_file(args, temp_dir, out_dir)
    model_path, ml_scale_path = check_model_files(args, temp_dir, out_dir) 
    names_of_columns = ['ssv', 'ssv1', 'ssv2', 'wsv', 'gtgt', 'gt', 'us', 'gs', 'aus', 'up_u', 'mid_u', 'tgta', 'r1-r2', 'r2-r3', 'rat_r1', 'rat_r2', 'rat_r3', 'no_r2_rat_r1', 'no_r2_rat_r3', 'small_r1-r2', 'small_r2-r3', 'small_rat_r1', 'small_rat_r2', 'small_rat_r3', 'small_no_r2_rat_r1', 'small_no_r2_rat_r3', 'stringtie']
    if machine_learning_mode:
        names_of_columns = names_of_columns+['lab']
    verbose = args.verbose
    if verbose:
        sys.stderr.write('Verbose enabled: Updates will be printed'+'\n')
        time_now('Starting aptardi analysis..', starting = True) 
    if verbose and debugging:
        sys.stderr.write('Debugging enabled: Intermediate files will be kept'+'\n')
    bed_list = stg_gtf_to_bed(stg_gtf, temp_dir)
    rec = format_bed_list(bed_list)
    bam_to_bedgraph(bam_file, temp_dir, bg_file)
    get_chrom_sizes(fasta_file, chr_sizes_name, temp_dir)
    chrom_dict = load_chrom_sizes(chr_sizes_name, temp_dir)
    extract_extended_and_orig_end_utr(rec, chrom_dict, temp_dir, 'ends.bed')
    extract_extension(rec, chrom_dict, temp_dir, 'extensions.bed')
    extract_start_utr(rec, temp_dir, 'starts.bed')
    intersect_extensions_start_utrs('extensions.bed', 'starts.bed', temp_dir, 'intersect_ext_start.bed')
    refine_utr_extension(temp_dir, 'intersect_ext_start.bed', 'ends.bed', 'refined.bed')
    utr_cov_dict, utr_seq_depth = extract_rna_seq_coverage_utr(bg_file, temp_dir, 'refined.bed')           
    end_base_dict = identify_transcript_end_base(utr_cov_dict)
    redefine_refined_utrs('refined.bed', end_base_dict, 'fin.bed', temp_dir)                                    
    window_refined_utrs(temp_dir, 'fin.bed', 'window.bed')              
    merge_refined_utrs(temp_dir, 'window.bed', 'merged.bed')
    define_data_instances(temp_dir, 'inst.bed', chrom_dict)
    tmp_fasta_file = check_fasta_headers(fasta_file, temp_dir, 'headers.txt')       
    get_sequence(fasta_file, 'inst.bed', 'seq.bed', temp_dir, tmp_fasta_file)
    inst_dict = load_data_instances('seq.bed', temp_dir)
    dna_feat_dict = extract_dna_features(inst_dict)     
    dna_freq_feat_dict = base_freq_features(inst_dict)       
    cov_dict, seq_depth = extract_rna_seq_coverage(bg_file, temp_dir, inst_dict)
    rna_feat_dict = extract_local_rna_seq_features(cov_dict, seq_depth)
    surr_rna_feat_dict = extract_rna_seq_features(cov_dict, seq_depth)                                      
    sg_ends_dict, sg_ends_minus_dict = stringtie_ends_as_dict(rec)
    sg_feat_dict = extract_stringtie_ends_as_feature(sg_ends_dict, sg_ends_minus_dict, inst_dict)
    label_plus_loc, label_minus_loc = extract_label_locations(polya_sites, temp_dir)
    lab_dict = extract_labels(label_plus_loc, label_minus_loc, inst_dict)
    ml_dict = combine_feat_lab_dicts(dna_feat_dict, dna_freq_feat_dict, surr_rna_feat_dict, rna_feat_dict, sg_feat_dict, lab_dict)
    ml_list = bi_di_format(ml_dict, names_of_columns)
    ml_dat, ml_names = split_names_data(ml_list)
    idx_con, idx_bi = get_continuous_vars(ml_dat)      
    if machine_learning_mode:
        X_train, X_test, y_train, y_test, train_idx, test_idx, feat_shape, ml_len, ml_scale = bi_di_format_features(ml_dat, idx_con, 99)
        weights = generate_weights(y_train, 99)
        class_weights = dist_class_weights(weights, y_train)
        model_cp = ModelCheckpoint(out_dir+mod_name, verbose = 0, save_best_only = True)
        model = Sequential()
        model.add(Masking(mask_value = 99, input_shape = (ml_len, feat_shape)))
        model.add(Bidirectional(LSTM(20, return_sequences = True), input_shape = (ml_len, feat_shape)))
        model.add(TimeDistributed(Dense(1, activation = 'sigmoid')))
        model.compile(loss = 'binary_crossentropy', optimizer = 'adam', metrics = [Precision(), Recall()])
        history = model.fit(X_train, y_train, validation_split = 0.25, epochs = 25, verbose = 0, class_weight = class_weights, callbacks = [model_cp])
        save_files([ml_scale], [sc_name], out_dir)
        new_define_data_instances(temp_dir, 'inst.bed', chrom_dict)
        get_sequence(fasta_file, 'inst.bed', 'seq.bed', temp_dir, tmp_fasta_file)
        inst_dict = load_data_instances('seq.bed', temp_dir)
        dna_feat_dict = extract_dna_features(inst_dict)   
        dna_freq_feat_dict = base_freq_features(inst_dict)       
        cov_dict, seq_depth = extract_rna_seq_coverage(bg_file, temp_dir, inst_dict)
        rna_feat_dict = extract_local_rna_seq_features(cov_dict, seq_depth)
        surr_rna_feat_dict = extract_rna_seq_features(cov_dict, seq_depth)         
        sg_ends_dict, sg_ends_minus_dict = stringtie_ends_as_dict(rec)
        sg_feat_dict = extract_stringtie_ends_as_feature(sg_ends_dict, sg_ends_minus_dict, inst_dict)
        ml_dict = combine_feat_lab_dicts(dna_feat_dict, dna_freq_feat_dict, surr_rna_feat_dict, rna_feat_dict, sg_feat_dict, lab_dict)
        ml_list = bi_di_format(ml_dict, names_of_columns)
        ml_dat, ml_names = split_names_data(ml_list)
    else:
        model = keras_load_mod(model_path)
        ml_scale = open_files(ml_scale_path)
        ml_len = model.layers[0].output_shape[1]
    X_test = format_features(ml_dat, idx_con, 99, ml_scale, ml_len)
    ml_names = trunc_names(ml_names, ml_len)
    ed = add_prob_pred(X_test, ml_len, model)
    preds = remove_mask_add_names(ed, 99, ml_names) 
    if isinstance(res_file, str):   
        in_fi = open(out_dir+res_file, 'w')         
        final_gtf(preds, stg_gtf, in_fi, threshold, temp_dir)
    else:
        final_gtf(preds, stg_gtf, res_file, threshold, temp_dir)
    try:
        os.rmdir(temp_dir)
    except OSError:
        pass
    try:
        os.rmdir(out_dir)
    except OSError:
        pass
    if verbose:
        time_now('', finishing = True) 
