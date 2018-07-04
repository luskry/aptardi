def chunkify(fname, size):
    
    fileEnd = os.path.getsize(fname)
    
    with open(fname, 'rb') as f:
        chunkEnd = f.tell()
        while True:
            chunkStart = chunkEnd
            f.seek(size, 1)
            f.readline()
            chunkEnd = f.tell()
            yield chunkStart, chunkEnd - chunkStart
            if chunkEnd > fileEnd:
                break
