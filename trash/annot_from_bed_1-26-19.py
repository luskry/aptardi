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
            elif strand != '.':
                pass   
            else:
                continue
            
            rel_exon_starts = (',').join([str(int(i)-int(tx_start)) for i in abs_exon_starts])
            exon_sizes = (',').join(str(i) for i in exon_sizes)
            write_line = [fields[0], tx_start, tx_end, fields[3], '0', strand, fields[6], fields[7], '255, 0, 0', fields[9], exon_sizes, rel_exon_starts]
            output_write.writelines('\t'.join(write_line) + '\n')
     
    output_write.close()
