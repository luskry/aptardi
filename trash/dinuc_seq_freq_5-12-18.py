def dinuc_seq_freq(UTR_bed_path, tmp_dir):
    
    UTR_dict = {}
    base_combos = ['AA', 'AC','AG','AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG','GT', 'TA', 'TC', 'TG', 'TT']
    upstream_region_one = 70
    upstream_region_two = 70
    downstream_region_one = 30
    downstream_region_two = 30
    
    with open(temp_dir+UTR_bed_path, 'r') as f:
        for line in f:
            fields = line.strip('\n').split('\t')
            UTR_dict[fields[3]] = fields[0:3]+[fields[5]]
            bases = fields[-1]
            for i in range(len(bases)):
                if i >= upstream_region_one+upstream_region_two and i < (len(bases) - (downstream_region_one+downstream_region_two)): 
                    upstream_region_one_probe = range(i - (upstream_region_one+upstream_region_two), i - (upstream_region_two))
                    upstream_region_two_probe = range(i - upstream_region_two, i)
                    downstream_region_one_probe = range(i, i + downstream_region_one)
                    downstream_region_two_probe = range(i + downstream_region_one, i + (downstream_region_one+downstream_region_two))
                    seq_list_upstream_region_one = []
                    seq_list_upstream_region_two = []
                    seq_list_downstream_region_one = []
                    seq_list_downstream_region_two = []
                    for nuc in upstream_region_one_probe:
                        first_base = bases[nuc]
                        second_base = bases[nuc + 1]
                        dinuc = first_base+second_base
                        seq_list_upstream_region_one.append(dinuc)
                        seq_list_upstream_region_one_final = [[x, seq_list_upstream_region_one.count(x)] for x in set(seq_list_upstream_region_one)]
                        for combo in base_combos:
                            if combo not in [x[0] for x in seq_list_upstream_region_one_final]:
                                seq_list_upstream_region_one_final.append([combo, 0])
                        seq_list_upstream_region_one_final = sorted(seq_list_upstream_region_one_final, key = lambda x: x[0])
                        dinuc_sum = sum(x[1] for x in seq_list_upstream_region_one_final)
                        dinuc_freq = [x[1]/dinuc_sum for x in seq_list_upstream_region_one_final]
                        seq_list_upstream_region_one_final = [x for x in seq_list_upstream_region_one_final if 'N' not in x[0]]
                        for sub in range(len(seq_list_upstream_region_one_final)):
                            seq_list_upstream_region_one_final[sub][1] = dinuc_freq[sub]
                    for nuc in upstream_region_two_probe:
                        first_base = bases[nuc]
                        second_base = bases[nuc + 1]
                        dinuc = first_base+second_base
                        seq_list_upstream_region_two.append(dinuc)
                        seq_list_upstream_region_two_final = [[x, seq_list_upstream_region_two.count(x)] for x in set(seq_list_upstream_region_two)]
                        for combo in base_combos:
                            if combo not in [x[0] for x in seq_list_upstream_region_two_final]:
                                seq_list_upstream_region_two_final.append([combo, 0])
                        seq_list_upstream_region_two_final = sorted(seq_list_upstream_region_two_final, key = lambda x: x[0])
                        dinuc_sum = sum(x[1] for x in seq_list_upstream_region_two_final)
                        dinuc_freq = [x[1]/dinuc_sum for x in seq_list_upstream_region_two_final]
                        seq_list_upstream_region_two_final = [x for x in seq_list_upstream_region_two_final if 'N' not in x[0]]
                        for sub in range(len(seq_list_upstream_region_two_final)):
                            seq_list_upstream_region_two_final[sub][1] = dinuc_freq[sub]
                    for nuc in downstream_region_one_probe:
                        first_base = bases[nuc]
                        second_base = bases[nuc + 1]
                        dinuc = first_base+second_base
                        seq_list_downstream_region_one.append(dinuc)
                        seq_list_downstream_region_one_final = [[x, seq_list_downstream_region_one.count(x)] for x in set(seq_list_downstream_region_one)]
                        for combo in base_combos:
                            if combo not in [x[0] for x in seq_list_downstream_region_one_final]:
                                seq_list_downstream_region_one_final.append([combo, 0])
                        seq_list_downstream_region_one_final = sorted(seq_list_downstream_region_one_final, key = lambda x: x[0])
                        dinuc_sum = sum(x[1] for x in seq_list_downstream_region_one_final)
                        dinuc_freq = [x[1]/dinuc_sum for x in seq_list_downstream_region_one_final]
                        seq_list_downstream_region_one_final = [x for x in seq_list_downstream_region_one_final if 'N' not in x[0]]
                        for sub in range(len(seq_list_downstream_region_one_final)):
                            seq_list_downstream_region_one_final[sub][1] = dinuc_freq[sub]
                    for nuc in downstream_region_two_probe:
                        first_base = bases[nuc]
                        second_base = bases[nuc + 1]
                        dinuc = first_base+second_base
                        seq_list_downstream_region_two.append(dinuc)
                        seq_list_downstream_region_two_final = [[x, seq_list_downstream_region_two.count(x)] for x in set(seq_list_downstream_region_two)]
                        for combo in base_combos:
                            if combo not in [x[0] for x in seq_list_downstream_region_two_final]:
                                seq_list_downstream_region_two_final.append([combo, 0])
                        seq_list_downstream_region_two_final = sorted(seq_list_downstream_region_two_final, key = lambda x: x[0])
                        dinuc_sum = sum(x[1] for x in seq_list_downstream_region_two_final)
                        dinuc_freq = [x[1]/dinuc_sum for x in seq_list_downstream_region_two_final]
                        seq_list_downstream_region_two_final = [x for x in seq_list_downstream_region_two_final if 'N' not in x[0]]
                        for sub in range(len(seq_list_downstream_region_two_final)):
                            seq_list_downstream_region_two_final[sub][1] = dinuc_freq[sub]
                    UTR_dict[fields[3]].append(seq_list_upstream_region_one_final+seq_list_upstream_region_two_final+seq_list_downstream_region_one_final+seq_list_downstream_region_two_final)
    return UTR_dict
