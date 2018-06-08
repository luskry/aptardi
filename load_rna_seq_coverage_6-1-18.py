def load_rna_seq_coverage(All_Wig_files, UTR_Annotation_file, tmp_dir):
    
    UTR_events_dict = {}
    All_Samples_Total_depth = []
    All_samples_extracted_3UTR_coverage_dict = {}
    All_samples_extracted_3UTR_coverage_dict_minus = {}
    upstream_region_one = 70
    upstream_region_two = 70
    downstream_region_one = 30
    downstream_region_two = 30
    region_chunks = 5
    
    with open(tmp_dir+UTR_Annotation_file,'r') as f: 
        for line in f:
            fields = line.strip('\n').split('\t')
            curr_chr = fields[0]
            if not curr_chr.startswith('chr'):
                curr_chr = 'chr' + curr_chr
            region_start = fields[1]
            region_end = fields[2]
            UTR_pos = '%s:%s-%s' % (curr_chr, region_start, region_end)
            UTR_events_dict[fields[3]] = [fields[0], region_start, region_end, fields[5], UTR_pos]
       
    for curr_wig_file in All_Wig_files:
        
        cur_sample_total_depth = 0
        num_line = 0
        curr_sample_All_chroms_coverage_dict = {}
        curr_sample_All_chroms_coverage_dict_minus = {}
        
        with open(tmp_dir+curr_wig_file, 'r') as f:
            for line in f:
                if '#' not in line:
                    fields = line.strip('\n').split('\t')
                    chrom_name = fields[0]
                    if not chrom_name.startswith('chr'):
                          chrom_name = 'chr' + chrom_name
                    region_start = int(fields[1])
                    region_end = int(fields[2])
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
        All_Samples_Total_depth.append(cur_sample_total_depth)
        curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)
        curr_sample_All_chroms_coverage_dict_minus[chrom_name][1].append(0)
        
        for curr_3UTR_event_id in UTR_events_dict:
            curr_3UTR_structure = UTR_events_dict[curr_3UTR_event_id]
            curr_chr = curr_3UTR_structure[0]
            
            if UTR_events_dict[curr_3UTR_event_id][-2] == '+':
                if curr_chr in curr_sample_All_chroms_coverage_dict:
                    curr_chr_coverage = curr_sample_All_chroms_coverage_dict[curr_chr]
                    region_start = int(curr_3UTR_structure[1])
                    region_end = int(curr_3UTR_structure[2])
                    left_region_index = bisect(curr_chr_coverage[0],region_start)
                    right_region_index = bisect(curr_chr_coverage[0],region_end)
                    extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]
                    extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                    extracted_3UTR_region.insert(0,region_start)
                    extracted_3UTR_region.append(region_end)
                    
                    if curr_3UTR_event_id not in All_samples_extracted_3UTR_coverage_dict:
                        All_samples_extracted_3UTR_coverage_dict[curr_3UTR_event_id] = []
                    All_samples_extracted_3UTR_coverage_dict[curr_3UTR_event_id].append([extracted_coverage,extracted_3UTR_region])
           
            elif UTR_events_dict[curr_3UTR_event_id][-2] == '-':
                if curr_chr in curr_sample_All_chroms_coverage_dict_minus:
                    curr_chr_coverage = curr_sample_All_chroms_coverage_dict_minus[curr_chr]
                    region_start = int(curr_3UTR_structure[1])
                    region_end = int(curr_3UTR_structure[2])
                    left_region_index = bisect(curr_chr_coverage[0],region_start)
                    right_region_index = bisect(curr_chr_coverage[0],region_end)
                    extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]
                    extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                    extracted_3UTR_region.insert(0,region_start)
                    extracted_3UTR_region.append(region_end)
                    
                    if curr_3UTR_event_id not in All_samples_extracted_3UTR_coverage_dict_minus:
                        All_samples_extracted_3UTR_coverage_dict_minus[curr_3UTR_event_id] = []
                    All_samples_extracted_3UTR_coverage_dict_minus[curr_3UTR_event_id].append([extracted_coverage,extracted_3UTR_region])
    
    try:
        os.remove(tmp_dir+UTR_Annotation_file)
        for curr_wig_file in All_Wig_files:
            os.remove(tmp_dir+curr_wig_file)
    except OSError:
        pass
    else:  
        UTR_events_dict_fin = {}
        for curr_3UTR_id in UTR_events_dict:
            curr_3UTR_all_samples_bp_coverage = []
            if UTR_events_dict[curr_3UTR_id][-2] == '-' and curr_3UTR_id in All_samples_extracted_3UTR_coverage_dict_minus:
                curr_3UTR_coverage_wig = All_samples_extracted_3UTR_coverage_dict_minus[curr_3UTR_id]
                for curr_sample_curr_3UTR_coverage_wig in curr_3UTR_coverage_wig: 
                    bp_coverage = np.zeros(curr_sample_curr_3UTR_coverage_wig[1][-1] - curr_sample_curr_3UTR_coverage_wig[1][0])
                    relative_start = curr_sample_curr_3UTR_coverage_wig[1][0]
                    for i in range(len(curr_sample_curr_3UTR_coverage_wig[0])):
                        curr_region_start = curr_sample_curr_3UTR_coverage_wig[1][i] - relative_start
                        curr_region_end = curr_sample_curr_3UTR_coverage_wig[1][i+1] - relative_start
                        bp_coverage[curr_region_start:curr_region_end] = curr_sample_curr_3UTR_coverage_wig[0][i]
                    bp_coverage = bp_coverage[::-1]
                    curr_3UTR_all_samples_bp_coverage.append(bp_coverage)
            elif UTR_events_dict[curr_3UTR_id][-2] == '+' and curr_3UTR_id in All_samples_extracted_3UTR_coverage_dict:
                curr_3UTR_coverage_wig = All_samples_extracted_3UTR_coverage_dict[curr_3UTR_id]
                for curr_sample_curr_3UTR_coverage_wig in curr_3UTR_coverage_wig: 
                    bp_coverage = np.zeros(curr_sample_curr_3UTR_coverage_wig[1][-1] - curr_sample_curr_3UTR_coverage_wig[1][0])
                    relative_start = curr_sample_curr_3UTR_coverage_wig[1][0]
                    for i in range(len(curr_sample_curr_3UTR_coverage_wig[0])):
                        curr_region_start = curr_sample_curr_3UTR_coverage_wig[1][i] - relative_start
                        curr_region_end = curr_sample_curr_3UTR_coverage_wig[1][i+1] - relative_start
                        bp_coverage[curr_region_start:curr_region_end] = curr_sample_curr_3UTR_coverage_wig[0][i]
                    curr_3UTR_all_samples_bp_coverage.append(bp_coverage)
            UTR_events_dict_fin[curr_3UTR_id] = curr_3UTR_all_samples_bp_coverage
    
    norm_UTR_events_dict_fin = {}
    for curr_event in UTR_events_dict_fin:
        for i in range(len(UTR_events_dict_fin[curr_event])):
            exp_list = []
            exper_wide_norm = []
            cov_array = UTR_events_dict_fin[curr_event][i]
            probe_region = cov_array[upstream_region_one+upstream_region_two:len(cov_array)-(downstream_region_one+downstream_region_two)]
            exper_wide_norm.append((sum(probe_region)/np.array(All_Samples_Total_depth)[i])/len(probe_region))
            for base in range(len(cov_array)):
                 if base >= upstream_region_one+upstream_region_two and base < (len(cov_array) - (downstream_region_one+downstream_region_two)): 
                     reg_one_exp = np.sum(cov_array[base-(upstream_region_one+upstream_region_two):base-upstream_region_two])/sum(probe_region)
                     up_reg_exp_probe = [int(upstream_region_two/region_chunks)] * region_chunks
                     error = int(upstream_region_two/region_chunks * region_chunks - sum(up_reg_exp_probe))
                     up_reg_exp_probe[-1] = up_reg_exp_probe[-1] + error
                     small_up_reg = []
                     for z in range(len(up_reg_exp_probe)):
                         if z != 5:
                             small_up_reg.append(np.sum(cov_array[base-upstream_region_one+sum(up_reg_exp_probe[:z]):base-upstream_region_one+sum(up_reg_exp_probe[:z+1])])/sum(probe_region))
                     reg_two_exp = np.sum(cov_array[base+downstream_region_one:base+downstream_region_one+downstream_region_two])/sum(probe_region)
                     down_reg_exp_probe = [int(downstream_region_one/region_chunks)] * region_chunks
                     error = int(downstream_region_one/region_chunks * region_chunks - sum(down_reg_exp_probe))
                     down_reg_exp_probe[-1] = down_reg_exp_probe[-1] + error
                     small_down_reg = []
                     for z in range(len(down_reg_exp_probe)+1):
                         if z != 0:
                             small_down_reg.append(np.sum(cov_array[base+sum(down_reg_exp_probe[:z-1]):base+sum(down_reg_exp_probe[:z])])/sum(probe_region))
                     exp_list.append([np.array(reg_one_exp), np.array(small_up_reg), np.array(small_down_reg), np.array(reg_two_exp)])
            norm_UTR_events_dict_fin[curr_event] = exp_list
            norm_UTR_events_dict_fin[curr_event].append(np.array(exper_wide_norm[i]/len(cov_array)))
    
    return norm_UTR_events_dict_fin
                                                     