def get_chrom_sizes(chrom_file, tmp_dir):
    
    chrom_sizes_dict = {}
    chrom_file = tmp_dir+chrom_file
    
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
