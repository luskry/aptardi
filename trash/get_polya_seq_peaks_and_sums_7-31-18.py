def get_polya_seq_peaks_and_sums(y, min_dist = 16, thres_abs = True, thres = 30, rel_thres_sum = 3):    

    if not thres_abs:
        thres = thres * (np.max(y) - np.min(y)) + np.min(y)
        
    min_dist = int(min_dist)
    thres_sum = thres * rel_thres_sum
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
                     & (np.greater(y, thres)))[0]

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
    
    peaks_fin = []
    
    for i in peaks:
        
        try:
            peaks_sum = np.sum(y[i-min_dist+1:i])+np.sum(y[i:i+min_dist])
        except IndexError:
            peaks_sum = np.sum(y[i-min_dist+1:i])+np.sum(y[i:])
        
        if peaks_sum >= thres_sum:
            
            peaks_fin.append([i, peaks_sum])
    
    return peaks_fin
