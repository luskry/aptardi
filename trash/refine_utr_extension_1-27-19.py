def refine_utr_extension(tmp_dir, intersected_extension_start_utr_file, extended_end_utr_file, refined_utr_file):
    
    extended_utr_subtract_overlap = {}
    overlap_extension = {}
    output_write = open(tmp_dir+refined_utr_file, 'w') 
    
    with open(tmp_dir+ extended_end_utr_file, 'r') as f:
        for line in f:
            if line.strip().split('\t')[0] == '\n' or line[0] == '#':
                continue
            else:   
                field = line.strip().split('\t')
                name = field[3]
                extended_utr_subtract_overlap[name] = field
   
    with open(tmp_dir+intersected_extension_start_utr_file, 'r') as f:
        for line in f:
            if line.strip().split('\t')[0] == '\n' or line[0] == '#':
                continue
            else:   
                field = line.strip().split('\t')
                name = field[3]
                strand = field[5]
                start = field[7]
                end = field[8]
                if name not in overlap_extension:
                    overlap_extension[name] = [field[0:6], strand, [start, end]]
                else:
                    overlap_extension[name][2].extend([start, end])      
                    
    for event in extended_utr_subtract_overlap:
        if event in overlap_extension:
            if overlap_extension[event][1] == '-':
                if int(max(overlap_extension[event][2])) < int(extended_utr_subtract_overlap[event][2]):
                    extended_utr_subtract_overlap[event][1] = max(overlap_extension[event][2])
                else:
                    extended_utr_subtract_overlap[event][1] = extended_utr_subtract_overlap[event][6]
            elif overlap_extension[event][1] == '+':
                if int(min(overlap_extension[event][2])) > int(extended_utr_subtract_overlap[event][1]):
                    extended_utr_subtract_overlap[event][2] = min(overlap_extension[event][2])
                else:
                     extended_utr_subtract_overlap[event][2] = extended_utr_subtract_overlap[event][7]
            extended_utr_subtract_overlap[event][10] = str(int(extended_utr_subtract_overlap[event][2]) - int(extended_utr_subtract_overlap[event][1])) 
            
        output_write.writelines('\t'.join(extended_utr_subtract_overlap[event]) + '\n')
    
    output_write.close()

    try:
        os.remove(tmp_dir+extended_end_utr_file)
    except OSError:
        pass
    try:
        os.remove(tmp_dir+intersected_extension_start_utr_file)
    except OSError:
        pass
    
