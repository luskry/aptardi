def load_polya_seq_coverage(All_Wig_files_polya, UTR_Annotation_file, tmp_dir):
    
    output_write = open(tmp_dir+'labels.txt', 'w')
    min_cov = 10
    UTR_events_dict = {}
    All_Samples_Total_depth = []
    All_samples_extracted_3UTR_coverage_dict = {}
    All_samples_extracted_3UTR_coverage_dict_minus = {}
    upstream_region_one = 70
    upstream_region_two = 70
    downstream_region_one = 30
    downstream_region_two = 30
    
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
       
    for curr_wig_file in All_Wig_files_polya:
        
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
        All_Samples_Total_depth.append(cur_sample_total_depth)
        curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)
        curr_sample_All_chroms_coverage_dict_minus[chrom_name][1].append(0)
        
        for curr_3UTR_event_id in UTR_events_dict:
            curr_3UTR_structure = UTR_events_dict[curr_3UTR_event_id]
            curr_chr = curr_3UTR_structure[0]
            
            if UTR_events_dict[curr_3UTR_event_id][-2] == '+':
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
           
            elif UTR_events_dict[curr_3UTR_event_id][-2] == '-':
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
    
    try:
        for curr_wig_file in All_Wig_files_polya:
            os.remove(tmp_dir+curr_wig_file)
    except OSError:
        pass
    else:  
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
            indicator_list = []
            for i in range(len(curr_3UTR_all_samples_bp_coverage)):
                cov_array = curr_3UTR_all_samples_bp_coverage[i]
                probe_region = cov_array[upstream_region_one+upstream_region_two:len(cov_array)-(downstream_region_one+downstream_region_two)]
                probe_region[probe_region < min_cov] = 0
                probe_region[probe_region >= min_cov] = 1
                indicator_list.append(' '.join([str(int(item)) for item in probe_region])) 
              
            output_write.writelines('\t'.join([str(curr_3UTR_id), ' '.join(indicator_list)])+'\n')
 
    output_write.close()
    
    if output_write.closed == False:
        pass
    else:           
        os.system('sort %s -o %s'  % (tmp_dir+'labels.txt', tmp_dir+'sorted_'+'labels.txt'))
        os.remove(tmp_dir+'labels.txt')
