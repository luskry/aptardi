def extract_start_utr(bed_list_raw, tmp_dir):
    
    output_write = open(tmp_dir+'start_utr.bed', 'w')
    
    for fields in bed_list_raw and fields[9] != 0:
        if '_' not in fields[0]:
            strand = fields[5]
            tx_start = fields[1]
            tx_end = fields[2]  
            exon_sizes = [int(i) for i in fields[10].strip(',').split(',')]         
            
            if strand == '+':
                tx_end = str(int(tx_start) + int(exon_sizes[0]))
                exon_sizes = str(exon_sizes[0])
            elif strand == '-':
                tx_start = str(int(tx_end) - int(exon_sizes[-1]))
                exon_sizes = str(exon_sizes[-1]) 
            else:
                continue
                
            write_line = [fields[0], tx_start, tx_end, fields[3], '0', strand]
            output_write.writelines('\t'.join(write_line) + '\n')
     
    output_write.close()
