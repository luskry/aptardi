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
