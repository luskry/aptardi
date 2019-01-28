def merge_refined_utrs(tmp_dir, refined_utr_file, merged_refined_utr_file):

    os.system('sort -k1,1 -k2,2n %s | bedtools merge -i - -s -c 4,5,6 -o collapse,distinct,distinct -delim "/" > %s' % (tmp_dir+refined_utr_file, tmp_dir+merged_refined_utr_file))
    
    try:
        os.remove(tmp_dir+refined_utr_file)
    except OSError:
        pass
