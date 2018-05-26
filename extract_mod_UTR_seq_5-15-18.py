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
