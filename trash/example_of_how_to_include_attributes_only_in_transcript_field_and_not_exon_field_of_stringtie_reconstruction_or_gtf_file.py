def stg_gtf_to_bed(gtf_file, tmp_dir):
    
    
    def bed_line(estart, eend, field, nline, fpkm, tpm):
        
        estp = estart[0] - 1
        eedp = eend[-1]
        geneid = re.findall(r'gene_name \"([\w\.]+)\"', field[8])
        if len(geneid) == 0:
            geneid = re.findall(r'ref_gene_name \"([\w\.]+)\"', field[8])
        if len(geneid) == 0:
            geneid = re.findall(r'gene_id \"([\w\.]+)\"', field[8])   
        if len(geneid) == 0:
            geneid = [field[8].split('gene_name ')[1].split()[0].replace(';', '')]
        if len(geneid) == 0:
            geneid = [field[8].split('gene_id ')[1].split()[0].replace(';', '')]
        if len(geneid) == 0:
            geneid = 'NA'
        transid = [field[8].split('transcript_id ')[1].split()[0].replace(';', '').replace('"', '')][0]
        transid2 = re.findall(r'reference_id \"([\w\.]+)\"', field[8])
        if len(transid2) == 0:
            transid2 = 'NA'
        else:
            transid2 = transid2[0]
        if not field[0].startswith('chr'):
            field[0] = 'chr' + field[0]
        line = [field[0], str(estp), str(eedp), transid, transid2, field[6], fpkm, tpm, str(len(estart))]
        seglen = [eend[i] - estart[i] + 1 for i in range(len(estart))]
        segstart = [estart[i] - estart[0] for i in range(len(estart))]
        strl = str(seglen[0])
        for i in range(1, len(seglen)):
            strl += ',' + str(seglen[i])
        strs = str(segstart[0])
        for i in range(1, len(segstart)):
            strs +=',' + str(segstart[i])
        line.extend([strl, strs])
        
        return line
    
    
    extracted_bed_line = []   
    estart = []
    eend = []
    fpkm = ''
    tpm = ''
    nline = 0
    prevfield = []
    prevtransid = ''
    with open(tmp_dir+gtf_file, 'r') as f:
        for lines in f:
            if lines.strip().split('\t')[0] == '\n' or lines[0] == '#':
                continue
            else:
                field = lines.strip().split('\t')
                field[8:] = [(' ').join(field[8:])]
                nline = nline + 1
                if nline == 1 and field[1] not in ['Cufflinks', 'StringTie']:
                    print('Warning: The GTF file may not have been produced by Cufflinks or StringTie and therefore may give erroneous results')
                    pass
                if field[2] != 'exon' and field[2] != 'transcript':
                    continue        
                transid = [field[8].split('transcript_id ')[1].split()[0].replace(';', '').replace('"', '')]
                if len(transid) > 0:
                    transid = transid[0]
                else:
                    transid = ''
                if field[2] == 'transcript' or (prevtransid != '' and transid != '' and transid != prevtransid):
                    if len(estart) != 0:
                        extracted_bed_line.append(bed_line(estart, eend, prevfield, nline, fpkm, tpm))
                    estart = []
                    eend = []
                    fpkm = ''
                    tpm = ''
                prevfield = field
                prevtransid = transid
                if field[2] == 'exon':
                    est = int(field[3])
                    eed = int(field[4])
                    estart += [est]
                    eend += [eed]
                if field[2] == 'transcript':
                    fpkmval = re.findall(r'FPKM \"([\w\.]+)\"', field[8])
                    tpmval = re.findall(r'TPM \"([\w\.]+)\"', field[8])
                if len(fpkmval) == 0:
                    fpkm = '0'
                else:
                    fpkm = fpkmval[0]
                if len(tpmval) == 0:
                    tpm = '0'
                else:
                    tpm = tpmval[0]
        if len(estart) != 0:
            extracted_bed_line.append(bed_line(estart, eend, prevfield, nline, fpkm, tpm))
                
    return extracted_bed_line
