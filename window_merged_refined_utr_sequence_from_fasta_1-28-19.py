def sequence_from_fasta(local_fasta_file, tmp_dir, window_merged_refined_utr_file, window_merged_refined_utr_file_with_sequence):
    
    if local_fasta_file.endswith('.gz'):
        os.system('gunzip -f %s' % (local_fasta_file))
        local_fasta_file = local_fasta_file[:-3]
    
    os.system('bedtools getfasta -fi %s -bed %s -name -bedOut -s > %s' % (tmp_dir+local_fasta_file, tmp_dir+window_merged_refined_utr_file, tmp_dir+window_merged_refined_utr_file_with_sequence))

    try:
        os.remove(tmp_dir+window_merged_refined_utr_file)
    except OSError:
        pass
