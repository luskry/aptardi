def combine_text_files(tmp_dir):
    
    dinuc = tmp_dir+'sorted_dinuc_seq_freq_features.txt'
    rna_seq = tmp_dir+'sorted_rna_seq_features.txt'
    ssv = tmp_dir+'sorted_strong_signal_variant_indicator.txt'
    wws = tmp_dir+'sorted_weak_signal_variant_indicator.txt'
    ca = tmp_dir+'sorted_ca_element_indicator.txt'
    labs = tmp_dir+'sorted_labels.txt'
    os.system('gzip -dc %s | join -t $"\t" - %s | join -t $"\t" - %s | join -t $"\t" - %s | join -t $"\t" - %s  | join -t $"\t" - %s | gzip -f -9 > %s' % (dinuc, rna_seq, ssv, wws, ca, labs, tmp_dir+'features_and_labels.txt'))
    for file in [dinuc, rna_seq, ssv, wws, ca, labs]:
        os.remove(file)
    