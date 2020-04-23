def format_gencode_gtf(gencode_gtf_file, tmp_dir, output_gencode_file):

    with open(tmp_dir+gencode_gtf_file, 'r') as f:
        
            line = []
            prevtransid = ''
            temp_line = []
            
            for lines in f:
                if lines.strip().split('\t')[0] == '\n' or lines[0] == '#' or (lines.strip().split('\t')[2] != 'transcript' and lines.strip().split('\t')[2] != 'exon'):
                    continue
                else:
                    original_field = lines.strip().split('\t')
                    rep_field = [(' ').join(original_field[8:])]
                    field = original_field[0:8]+rep_field
                    transid = field[8].split('transcript_id ')[1].split()[0].replace(';', '')
                    if field[6] == '-' and field[2] == 'exon':
                        temp_line.append(original_field)
                    if transid != prevtransid:
                        temp_line = temp_line[::-1]
                        for i in temp_line:
                            line.append(i)
                        temp_line = []
                    if field[6] == '-' and field[2] == 'transcript':
                        line.append(original_field)
                    if  field[6] == '+':
                        line.append(original_field)
                    prevtransid = transid
            
            if field[6] == '-':
                temp_line = temp_line[::-1]
                for i in temp_line:
                    line.append(i)
                          
    with open(tmp_dir+output_gencode_file, 'w') as f:
        f.writelines('\t'.join(i) + '\n' for i in line)       
        
