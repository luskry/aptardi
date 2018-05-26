def sequence_from_fasta(local_fasta_file, tmp_dir):
    
    if local_fasta_file.endswith('.gz'):
        os.system('gunzip -f %s' % (local_fasta_file))
        local_fasta_file = local_fasta_file[:-3]
    
    os.system('bedtools getfasta -fi %s -bed %s -name -bedOut -s > %s' % (local_fasta_file, tmp_dir+'modified_annot.bed', tmp_dir+'modified_bed_with_sequence.bed'))

    try:
        os.remove(tmp_dir+'modified_annot.bed')
    except OSError:
        pass
