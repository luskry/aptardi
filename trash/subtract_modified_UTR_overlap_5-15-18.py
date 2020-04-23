def subtract_modified_UTR_overlap(UTR_bed_path, tmp_dir):
    
    os.system('sort -k1,1 -k2,2n -k3,3n -k6,6 -u %s > %s' % (tmp_dir+UTR_bed_path, tmp_dir+'UTR_duplicate_subtracted.bed'))
    
    try:
        os.remove(tmp_dir+UTR_bed_path)
    except OSError:
        pass
