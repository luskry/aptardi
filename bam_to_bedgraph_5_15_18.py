def bam_to_bedgraph(rna_seq_bam_path, polya_seq_bam_path, tmp_dir):
      
    os.system('bedtools genomecov -bg -strand + -ibam %s > %s' % (rna_seq_bam_path, tmp_dir+'rna_seq_coverage_plus_strand.bedgraph'))
    os.system('bedtools genomecov -bg -strand - -ibam %s > %s' % (rna_seq_bam_path, tmp_dir+'rna_seq_coverage_minus_strand.bedgraph'))
    os.system('bedtools unionbedg -i %s %s > %s' % (tmp_dir+'rna_seq_coverage_plus_strand.bedgraph', tmp_dir+'rna_seq_coverage_minus_strand.bedgraph', tmp_dir+'gold_standard_rna_seq_coverage.bedgraph'))
        
    try:
        os.remove(tmp_dir+'rna_seq_coverage_plus_strand.bedgraph')
    except OSError:
        pass
    
    try:
        os.remove(tmp_dir+'rna_seq_coverage_minus_strand.bedgraph')
    except OSError:
        pass
      
    os.system('bedtools genomecov -bg -strand + -5 -ibam %s > %s' % (polya_seq_bam_path, tmp_dir+'polya_seq_coverage_plus_strand.bedgraph'))
    os.system('bedtools genomecov -bg -strand - -5 -ibam %s > %s' % (polya_seq_bam_path, tmp_dir+'polya_seq_coverage_minus_strand.bedgraph'))
    os.system('bedtools unionbedg -i %s %s > %s' % (tmp_dir+'polya_seq_coverage_plus_strand.bedgraph', tmp_dir+'polya_seq_coverage_minus_strand.bedgraph', tmp_dir+'gold_standard_polya_seq_coverage.bedgraph'))
    
    try:
        os.remove(tmp_dir+'polya_seq_coverage_plus_strand.bedgraph')
    except OSError:
        pass
    
    try:
        os.remove(tmp_dir+'polya_seq_coverage_minus_strand.bedgraph')
    except OSError:
        pass
     