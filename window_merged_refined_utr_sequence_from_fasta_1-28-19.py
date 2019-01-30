def sequence_from_fasta(local_fasta_file, tmp_dir, window_merged_refined_utr_file, window_merged_refined_utr_file_with_sequence):
    
    window_size = 100
    signal_window_size = 36
    
    if local_fasta_file.endswith('.gz'):
        os.system('gunzip -f %s' % (tmp_dir+local_fasta_file))
        local_fasta_file = local_fasta_file[:-3]
    
    if window_size > signal_window_size:
        os.system('bedtools getfasta -fi %s -bed %s -name -bedOut -s > %s' % (tmp_dir+local_fasta_file, tmp_dir+window_merged_refined_utr_file, tmp_dir+window_merged_refined_utr_file_with_sequence))
    elif window_size < signal_window_size:
        get_sequence = signal_window_size - window_size
        os.system('awk \'BEGIN{FS=OFS="\t"} $6=="+" {$2 = $2 - %s;}1\' %s |  awk \'BEGIN{FS=OFS="\t"} $6=="-" {$3 = $3 + %s;}1\' | bedtools getfasta -fi %s -bed - -name -bedOut -s > %s' % (get_sequence, tmp_dir+window_merged_refined_utr_file, get_sequence, tmp_dir+local_fasta_file, tmp_dir+window_merged_refined_utr_file_with_sequence))

    try:
        os.remove(tmp_dir+window_merged_refined_utr_file)
    except OSError:
        pass
