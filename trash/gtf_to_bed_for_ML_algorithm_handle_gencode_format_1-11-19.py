def gtf_to_bed_for_ML_algorithm(gtf_file, tmp_dir):
    
    
    def bed_line(estart, eend, field, nline):
        
        allids = {}
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
        else:
            geneid = geneid[0]
        transid = re.findall(r'reference_id \"([\w\.]+)\"', field[8])
        if len(transid) == 0:
            transid = re.findall(r'transcript_id \"([\w\.]+)\"', field[8])
        if len(transid) == 0:
            transid = [field[8].split('transcript_id ')[1].split()[0].replace(';', '')]
        if len(transid) == 0:
            transid = 'Trans_' + str(nline)
        else:
            transid = transid[0]
        if transid in allids.keys():
            transid2 = transid + '_DUP' + str(allids[transid])
            allids[transid] = allids[transid] + 1
            transid = transid2
        else:
            allids[transid] = 1
        if not field[0].startswith('chr'):
            field[0] = 'chr' + field[0]
        line = [field[0], str(estp), str(eedp), '|'.join([transid, geneid, field[0], str(estp), str(eedp), field[6]]), '0', field[6], str(estp), str(eedp), '255,0,0', str(len(estart))]
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
    nline = 0
    prevfield = []
    prevtransid = ''
    with open(gtf_file, 'r') as f:
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
                transid = re.findall(r'transcript_id \"([\w\.]+)\"', field[8])
                if len(transid) == 0:
                    transid = [field[8].split('transcript_id ')[1].split()[0].replace(';', '')]
                if len(transid) > 0:
                    transid = transid[0]
                else:
                    transid = ''
                if field[2] == 'transcript' or (prevtransid != '' and transid != '' and transid != prevtransid):
                    if len(estart) != 0:
                        extracted_bed_line.append(bed_line(estart, eend, prevfield, nline))
                    estart = []
                    eend = []
                prevfield = field
                prevtransid = transid
                if field[2] == 'exon':
                    est = int(field[3])
                    eed = int(field[4])
                    estart += [est]
                    eend += [eed]
        if len(estart) != 0:
            extracted_bed_line.append(bed_line(estart, eend, prevfield, nline))
                
    return extracted_bed_line

    if os.path.dirname(gtf_file) == tmp_dir[:-1]:
    
        try:
            os.remove(gtf_file)
        except OSError:
            pass
        
