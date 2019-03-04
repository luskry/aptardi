def quote_gtf(input_gtf, output_gtf, tmp_dir):

    output_write = open(tmp_dir+output_gtf, 'w')
        
    with open(tmp_dir+input_gtf, 'r') as f:
        for line in f:
            if line.strip().split('\t')[0] == '\n' or line[0] == '#':
                continue
            else:
                fields = line.strip('\n').split('\t')
                key = re.findall('\w+(?=\s* )', fields[8])
                att = re.findall('[\w.\-\_]+(?=\s*;)', fields[8])
                att_rep = []
                for i in att:
                    att_rep.append('"%s"' % (i))
                fields[8] = ('').join([m+' '+n+'; ' for m, n in zip(key, att_rep)]).strip()
        
            output_write.writelines('\t'.join(fields) + '\n')
     
    output_write.close()
