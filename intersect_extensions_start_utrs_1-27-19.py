def intersect_extensions_start_utrs(extension_file, start_utr_file, tmp_dir, intersected_extension_start_utr_file):
    
    os.system('bedtools intersect -a %s -b %s -s -wo > %s' % (tmp_dir+extension_file, temp_dir+start_utr_file, tmp_dir+intersected_extension_start_utr_file))

    try:
        os.remove(tmp_dir+extension_file)
    except OSError:
        pass
    try:
        os.remove(temp_dir+start_utr_file)
    except OSError:
        pass
         