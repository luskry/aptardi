def window_merged_refined_utrs(tmp_dir, merged_refined_utr_file, window_merged_refined_utr_file):
    
    window_size = 100
    output_write = open(tmp_dir+window_merged_refined_utr_file, 'w') 
    
    with open(tmp_dir+ merged_refined_utr_file, 'r') as f:
        for line in f:
            if line.strip().split('\t')[0] == '\n' or line[0] == '#' or (int(line.strip().split('\t')[2])-int(line.strip().split('\t')[1]))/window_size < 3:
                continue
            else:
                field = line.strip().split('\t')
                start = int(field[1])
                end = int(field[2])
                utr_size = str(end-start-(end-start)%window_size)
                if field[5] == '+':
                    field[2] = str(end-(end-start)%window_size)
                elif field[5] == '-':
                    field[1] = str(start+(end-start)%window_size)
                field.append(utr_size)
            
                output_write.writelines('\t'.join(field) + '\n')
        output_write.close()
        
    try:
        os.remove(tmp_dir+merged_refined_utr_file)
    except OSError:
        pass
      