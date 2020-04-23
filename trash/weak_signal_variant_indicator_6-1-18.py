def weak_signal_variant_indicator(UTR_bed_path, tmp_dir):
    
    UTR_dict = {}
    signal_variants = ['AAGAAA', 'AAAAAG', 'AATACA', 'TATAAA', 'ACTAAA', 'AGTAAA', 'GATAAA', 'AATATA', 'CATAAA', 'AATAGA']
    upstream_region_one = 70
    upstream_region_two = 70
    downstream_region_one = 30
    downstream_region_two = 30
    
    with open(tmp_dir+UTR_bed_path, 'r') as f:
        for line in f:
            fields = line.strip('\n').split('\t')
            bases = fields[-1]
            i_line = -1
            signal_variant_indicator = np.zeros(len(bases) - (upstream_region_one + upstream_region_two + downstream_region_one + downstream_region_two))
            for i in range(len(bases)):
                if i >= (upstream_region_one+upstream_region_two) and i < (len(bases) - (downstream_region_one+downstream_region_two)): 
                    probe = range(i - 30, i - 5) #window = 36th-7th base upstream base of interest not including base of interest
                    hex_seq_list = []
                    i_line += 1
                    for nuc in probe:
                        hexamer = bases[nuc-6:nuc]
                        hex_seq_list.append(hexamer)
                        if any(i in signal_variants for i in hex_seq_list):
                            signal_variant_indicator[i_line] = 1
                        else:
                            signal_variant_indicator[i_line] = 0
            UTR_dict[fields[3]] = signal_variant_indicator
            
    return UTR_dict
