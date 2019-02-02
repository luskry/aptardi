def extract_label_locations(gencode_manually_annotated_file, tmp_dir):
    
    label_dict = {}
    label_dict_minus = {}
    
    with open(tmp_dir+gencode_manually_annotated_file, 'r') as f:
            for line in f:
                if line.strip().split('\t')[0] == '\n' or line[0] == '#':
                    continue
                else:   
                    fields = line.strip().split('\t')
                    chrom_name = fields[0]
                    strand = fields[3]
                    polya_site = int(fields[1])
                    if not chrom_name.startswith('chr'):
                        chrom_name = 'chr' + chrom_name
                    if strand == '-':
                        polya_site = polya_site - 1
                        if chrom_name not in label_dict_minus:
                            label_dict_minus[chrom_name] = [polya_site]
                        else:
                            label_dict_minus[chrom_name].extend([polya_site])
                    elif strand == '+':
                        if chrom_name not in label_dict:
                            label_dict[chrom_name] = [polya_site]
                        else:
                            label_dict[chrom_name].extend([polya_site])
    
    return label_dict, label_dict_minus
