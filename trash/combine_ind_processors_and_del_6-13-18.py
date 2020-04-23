def combine_ind_processors_and_del(ind_proc_file_list, tmp_dir, combined_file):
    
    with gzip.open(tmp_dir+combined_file, 'wt') as outfile:
        for fname in ind_proc_file_list:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
            if infile.closed == False:
                pass
            else:
                os.remove(fname)
    os.system('gzip -dc %s | sort | gzip -f -9 > %s'  % (tmp_dir+combined_file, tmp_dir+'sorted_'+combined_file))         
    os.remove(tmp_dir+combined_file)
 
    
combine_ind_processors_and_del(files_list_dinuc_seq_freq, temp_dir, 'dinuc_seq_freq.txt')
