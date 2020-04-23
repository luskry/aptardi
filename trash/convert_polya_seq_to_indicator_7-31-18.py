def convert_polya_seq_to_indicator(y, min_dist = 16):    
        
    min_dist = int(min_dist)
    dy = np.diff(y)
    zeros,=np.where(dy == 0)
    
    if len(zeros) == len(y) - 1:
        return np.array([])
    
    while len(zeros):
      
        zerosr = np.hstack([dy[1:], 0.])
        zerosl = np.hstack([0., dy[:-1]])
        dy[zeros]=zerosr[zeros]
        zeros,=np.where(dy == 0)
        dy[zeros]=zerosl[zeros]
        zeros,=np.where(dy == 0)

    peaks = np.where((np.hstack([dy, 0.]) < 0.)
                     & (np.hstack([0., dy]) > 0.)
                     & (np.greater(y, 0)))[0]

    if peaks.size > 1 and min_dist > 1:
        highest = peaks[np.argsort(y[peaks])][::-1]
        rem = np.ones(y.size, dtype = bool)
        rem[peaks] = False

        for peak in highest:
            if not rem[peak]:
                sl = slice(max(0, peak - min_dist), peak + min_dist + 1)
                rem[sl] = True
                rem[peak] = False

        peaks = np.arange(y.size)[~rem]
    
    indicator_fin = np.zeros(len(y))
    np.put(indicator_fin, peaks, 1)
        
    return indicator_fin
