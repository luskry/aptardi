def ca_element_indicator(UTR_bed_path, tmp_dir):
    
    output_write = open(tmp_dir+'ca_element_indicator.txt', 'w')
    upstream_region_one = 70
    upstream_region_two = 70
    downstream_region_one = 30
    downstream_region_two = 30
    
    with open(tmp_dir+UTR_bed_path, 'r') as f:
        for line in f:
            fields = line.strip('\n').split('\t')
            bases = fields[-1]
            i_line = -1
            ca_element_indicator = np.zeros(len(bases)-(upstream_region_one+upstream_region_two+downstream_region_one+downstream_region_two))
            for i in range(len(bases)):
                if i >= (upstream_region_one+upstream_region_two) and i < (len(bases)-(downstream_region_one+downstream_region_two)): 
                    i_line += 1
                    first_base = bases[i - 2]
                    second_base = bases[i - 1]
                    dinuc = first_base+second_base
                    if dinuc == 'CA':
                        ca_element_indicator[i_line] = 1
                    else:
                        ca_element_indicator[i_line] = 0
            output_write.writelines('\t'.join([str(fields[3]), ' '.join([str(int(item)) for item in ca_element_indicator])])+'\n')

    output_write.close()
    
    if output_write.closed == False:
        pass
    else:           
        os.system('sort %s -o %s'  % (tmp_dir+'ca_element_indicator.txt', tmp_dir+'sorted_'+'ca_element_indicator.txt'))
        os.remove(tmp_dir+'ca_element_indicator.txt')
                 