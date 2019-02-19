def extract_features_and_labels(seq_window_merged_refined_utrs_file, tmp_dir, bedgraph_files, gencode_label_dict, gencode_label_dict_minus):
    
    window_size = 100
    signal_window_size = 36
    min_upstream_signal_size = 7
    UTR_events_dict = {}
    All_Samples_Total_depth = []
    All_samples_extracted_3UTR_coverage_dict = {}
    All_samples_extracted_3UTR_coverage_dict_minus = {}
    features_labels_dict = {}

    with open(tmp_dir+seq_window_merged_refined_utrs_file, 'r') as f:
            for line in f:
                fields = line.strip('\n').split('\t')
                name = fields[3]
                start = int(fields[1])
                end = int(fields[2])
                bases = fields[-1]
                
                for idx, i in enumerate(range(start, end - window_size * 2, window_size)):
                    region_start = i
                    region_end = i + window_size * 3
                    rel_start = idx * window_size
                    rel_end = (idx + 3) * window_size
                    utr_pos = '%s-%s' % (i + window_size, i + window_size * 2)
                    curr_bases = bases[rel_start + window_size - signal_window_size:rel_end - window_size - min_upstream_signal_size]
                    UTR_events_dict['%s|%s-%s' % (name, i + window_size, i + window_size * 2)] = [fields[0], region_start, region_end, fields[5], utr_pos, curr_bases]      
    
    for curr_3UTR_event_id in UTR_events_dict:
        features_labels_dict[curr_3UTR_event_id] = []
        curr_3UTR_structure = UTR_events_dict[curr_3UTR_event_id]
        curr_chr = curr_3UTR_structure[0]
        bases = UTR_events_dict[curr_3UTR_event_id][-1]
        features_labels_dict[curr_3UTR_event_id].append(signal_variant_indicator(bases, ['AATAAA', 'ATTAAA']))
        features_labels_dict[curr_3UTR_event_id].append(signal_variant_indicator(bases, ['AAGAAA', 'AAAAAG', 'AATACA', 'TATAAA', 'ACTAAA', 'AGTAAA', 'GATAAA', 'AATATA', 'CATAAA', 'AATAGA']))                  
    
        if UTR_events_dict[curr_3UTR_event_id][-3] == '+':              
            if curr_chr in gencode_label_dict:
                probe_start, probe_end = int(UTR_events_dict[curr_3UTR_event_id][-2].split('-')[0]), int(UTR_events_dict[curr_3UTR_event_id][-2].split('-')[1])
                if any([i for i in gencode_label_dict[curr_chr] if i >= probe_start and i < probe_end]):
                    features_labels_dict[curr_3UTR_event_id].append(1)
                else:
                    features_labels_dict[curr_3UTR_event_id].append(0)
            else:
                features_labels_dict[curr_3UTR_event_id].append(0)   
        elif UTR_events_dict[curr_3UTR_event_id][-3] == '-':              
            if curr_chr in gencode_label_dict_minus:
                probe_start, probe_end = int(UTR_events_dict[curr_3UTR_event_id][-2].split('-')[0]), int(UTR_events_dict[curr_3UTR_event_id][-2].split('-')[1])
                if any([i for i in gencode_label_dict_minus[curr_chr] if i >= probe_start and i < probe_end]):
                    features_labels_dict[curr_3UTR_event_id].append(1)
                else:
                    features_labels_dict[curr_3UTR_event_id].append(0)
            else:
                features_labels_dict[curr_3UTR_event_id].append(0)           
    
    for curr_bedgraph in bedgraph_files:
        cur_sample_total_depth = 0
        num_line = 0
        curr_sample_All_chroms_coverage_dict = {}
        curr_sample_All_chroms_coverage_dict_minus = {}
        curr_chr = curr_3UTR_structure[0]
        
        with open(tmp_dir+curr_bedgraph, 'r') as f:
            for line in f:
                if '#' not in line:
                    fields = line.strip('\n').split('\t')
                    chrom_name = fields[0]
                    if not chrom_name.startswith('chr'):
                          chrom_name = 'chr' + chrom_name
                    region_start = int(float(fields[1]))
                    region_end = int(float(fields[2]))
                    cur_sample_total_depth += int(float(fields[-2])) * (region_end - region_start) + int(float(fields[-1])) * (region_end - region_start)
                    
                    if chrom_name not in curr_sample_All_chroms_coverage_dict:
                        curr_sample_All_chroms_coverage_dict[chrom_name] = [[0],[0]]
                    if region_start > curr_sample_All_chroms_coverage_dict[chrom_name][0][-1]:
                        curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_start)
                        curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)
                    curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_end)
                    curr_sample_All_chroms_coverage_dict[chrom_name][1].append(int(float(fields[-2])))
            
                    if chrom_name not in curr_sample_All_chroms_coverage_dict_minus:
                        curr_sample_All_chroms_coverage_dict_minus[chrom_name] = [[0],[0]]
                    if region_start > curr_sample_All_chroms_coverage_dict_minus[chrom_name][0][-1]:
                        curr_sample_All_chroms_coverage_dict_minus[chrom_name][0].append(region_start)
                        curr_sample_All_chroms_coverage_dict_minus[chrom_name][1].append(0)
                    curr_sample_All_chroms_coverage_dict_minus[chrom_name][0].append(region_end)
                    curr_sample_All_chroms_coverage_dict_minus[chrom_name][1].append(int(float(fields[-1])))
                num_line += 1
        if cur_sample_total_depth != 0:
            All_Samples_Total_depth.append(cur_sample_total_depth)
        else:
            print('Error: No read coverage detected for %s' % (curr_bedgraph))
            continue
        curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)
        curr_sample_All_chroms_coverage_dict_minus[chrom_name][1].append(0)    
    
        for curr_3UTR_event_id in UTR_events_dict:
            curr_3UTR_structure = UTR_events_dict[curr_3UTR_event_id]
            curr_chr = curr_3UTR_structure[0]
            curr_3UTR_all_samples_bp_coverage = []
            exp_list = []               
            
            if UTR_events_dict[curr_3UTR_event_id][-3] == '+':
                if curr_chr in curr_sample_All_chroms_coverage_dict:
                    curr_chr_coverage = curr_sample_All_chroms_coverage_dict[curr_chr]
                    region_start = int(float(curr_3UTR_structure[1]))
                    region_end = int(float(curr_3UTR_structure[2]))
                    left_region_index = bisect(curr_chr_coverage[0],region_start)
                    right_region_index = bisect(curr_chr_coverage[0],region_end)
                    extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]
                    extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                    extracted_3UTR_region.insert(0,region_start)
                    extracted_3UTR_region.append(region_end)
                    
                    if curr_3UTR_event_id not in All_samples_extracted_3UTR_coverage_dict:
                        All_samples_extracted_3UTR_coverage_dict[curr_3UTR_event_id] = []
                    All_samples_extracted_3UTR_coverage_dict[curr_3UTR_event_id].append([extracted_coverage,extracted_3UTR_region])
           
            elif UTR_events_dict[curr_3UTR_event_id][-3] == '-':
                if curr_chr in curr_sample_All_chroms_coverage_dict_minus:
                    curr_chr_coverage = curr_sample_All_chroms_coverage_dict_minus[curr_chr]
                    region_start = int(float(curr_3UTR_structure[1]))
                    region_end = int(float(curr_3UTR_structure[2]))
                    left_region_index = bisect(curr_chr_coverage[0],region_start)
                    right_region_index = bisect(curr_chr_coverage[0],region_end)
                    extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]
                    extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                    extracted_3UTR_region.insert(0,region_start)
                    extracted_3UTR_region.append(region_end)
                    
                    if curr_3UTR_event_id not in All_samples_extracted_3UTR_coverage_dict_minus:
                        All_samples_extracted_3UTR_coverage_dict_minus[curr_3UTR_event_id] = []
                    All_samples_extracted_3UTR_coverage_dict_minus[curr_3UTR_event_id].append([extracted_coverage,extracted_3UTR_region])
    
            if UTR_events_dict[curr_3UTR_event_id][-3] == '-' and curr_3UTR_event_id in All_samples_extracted_3UTR_coverage_dict_minus:
                curr_3UTR_coverage_wig = All_samples_extracted_3UTR_coverage_dict_minus[curr_3UTR_event_id]
                for curr_sample_curr_3UTR_coverage_wig in curr_3UTR_coverage_wig: 
                    bp_coverage = np.zeros(curr_sample_curr_3UTR_coverage_wig[1][-1] - curr_sample_curr_3UTR_coverage_wig[1][0])
                    relative_start = curr_sample_curr_3UTR_coverage_wig[1][0]
                    for i in range(len(curr_sample_curr_3UTR_coverage_wig[0])):
                        curr_region_start = curr_sample_curr_3UTR_coverage_wig[1][i] - relative_start
                        curr_region_end = curr_sample_curr_3UTR_coverage_wig[1][i+1] - relative_start
                        bp_coverage[curr_region_start:curr_region_end] = curr_sample_curr_3UTR_coverage_wig[0][i]
                    bp_coverage = bp_coverage[::-1]
                    curr_3UTR_all_samples_bp_coverage.append(bp_coverage)
            elif UTR_events_dict[curr_3UTR_event_id][-3] == '+' and curr_3UTR_event_id in All_samples_extracted_3UTR_coverage_dict:
                curr_3UTR_coverage_wig = All_samples_extracted_3UTR_coverage_dict[curr_3UTR_event_id]
                for curr_sample_curr_3UTR_coverage_wig in curr_3UTR_coverage_wig: 
                    bp_coverage = np.zeros(curr_sample_curr_3UTR_coverage_wig[1][-1] - curr_sample_curr_3UTR_coverage_wig[1][0])
                    relative_start = curr_sample_curr_3UTR_coverage_wig[1][0]
                    for i in range(len(curr_sample_curr_3UTR_coverage_wig[0])):
                        curr_region_start = curr_sample_curr_3UTR_coverage_wig[1][i] - relative_start
                        curr_region_end = curr_sample_curr_3UTR_coverage_wig[1][i+1] - relative_start
                        bp_coverage[curr_region_start:curr_region_end] = curr_sample_curr_3UTR_coverage_wig[0][i]
                    curr_3UTR_all_samples_bp_coverage.append(bp_coverage)              
             
            for i in range(len(curr_3UTR_all_samples_bp_coverage)):
                cov_array = curr_3UTR_all_samples_bp_coverage[i]
                
                reg_one_sum = (np.sum(cov_array[0:window_size])*1000000)/(np.array(All_Samples_Total_depth)[i]/1000000)
                reg_two_sum = (np.sum(cov_array[window_size:2*window_size])*1000000)/(np.array(All_Samples_Total_depth)[i]/1000000)
                reg_three_sum = (np.sum(cov_array[2*window_size:3*window_size])*1000000)/(np.array(All_Samples_Total_depth)[i]/1000000)
                
                reg_one_mean = (np.mean(cov_array[0:window_size])*1000000)/(np.array(All_Samples_Total_depth)[i]/1000000)
                reg_two_mean = (np.mean(cov_array[window_size:2*window_size])*1000000)/(np.array(All_Samples_Total_depth)[i]/1000000)
                reg_three_mean = (np.mean(cov_array[2*window_size:3*window_size])*1000000)/(np.array(All_Samples_Total_depth)[i]/1000000)
                
                reg_one_med = (np.median(cov_array[0:window_size])*1000000)/(np.array(All_Samples_Total_depth)[i]/1000000)
                reg_two_med = (np.median(cov_array[window_size:2*window_size])*1000000)/(np.array(All_Samples_Total_depth)[i]/1000000)
                reg_three_med = (np.median(cov_array[2*window_size:3*window_size])*1000000)/(np.array(All_Samples_Total_depth)[i]/1000000)
                
                un_norm_two_one_diff_med =  abs(reg_two_med-reg_one_med)
                un_norm_two_one_diff_mean =  abs(reg_two_mean-reg_one_mean)
                un_norm_three_two_diff_med = abs(reg_two_med-reg_three_med)
                un_norm_three_two_diff_mean = abs(reg_two_mean-reg_three_mean)
                
                if reg_one_sum+reg_two_sum != 0:
                    two_one_diff_med =  abs(reg_two_med-reg_one_med)/(reg_one_sum+reg_two_sum)
                    two_one_diff_mean =  abs(reg_two_mean-reg_one_mean)/(reg_one_sum+reg_two_sum)   
                else:
                    two_one_diff_med = 0
                    two_one_diff_mean = 0
                if reg_two_sum+reg_three_sum != 0:
                    three_two_diff_med =  abs(reg_two_med-reg_three_med)/(reg_three_sum+reg_two_sum)
                    three_two_diff_mean =  abs(reg_two_mean-reg_three_mean)/(reg_three_sum+reg_two_sum)
                    
                else:
                    three_two_diff_med = 0
                    three_two_diff_mean = 0
                    
                overall_diff_med = two_one_diff_med+three_two_diff_med
                overall_diff_mean = two_one_diff_mean+three_two_diff_mean
                un_norm_overall_diff_med = un_norm_two_one_diff_med+un_norm_three_two_diff_med
                un_norm_overall_diff_mean =un_norm_two_one_diff_mean+un_norm_three_two_diff_mean
                
                exp_list.append(overall_diff_med)
                
            features_labels_dict[curr_3UTR_event_id].insert(-1, sum(exp_list)) 
    
    return features_labels_dict
