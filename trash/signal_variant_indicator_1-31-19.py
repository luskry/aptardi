def signal_variant_indicator(bases, signal_variants):
    
    bases = bases.upper()
    hex_seq_list = []
    for i in range(len(bases) - 5):
        hex_seq_list.append(bases[i:i + 6])
    if any(sig in signal_variants for sig in hex_seq_list):
        return 1
    else:
        return 0
