def combine_feature_dicts(*arg):

    comb_dicts = {k:[d[k] for d in arg] for k in arg[0]}
    
    region_chunks = 5
    arr_dim = region_chunks*2+2
    fin_arr_all = np.array([]).reshape(0, 80)
    
    for event in comb_dicts:
        fin_arr = np.array([]).reshape(0, arr_dim)
        for features in range(len(comb_dicts[event])):
            if features == 0:
                for feature in range(len(comb_dicts[event][features])):
                    if feature != 1800:
                        feat_arr = np.hstack([x for x in comb_dicts[event][features][feature]])                        
                        fin_arr = np.vstack((fin_arr, feat_arr))
                    else:
                        overall_exp = np.repeat(comb_dicts[event][features][feature], len(comb_dicts[event][features])-1)
                fin_arr = np.hstack((fin_arr, overall_exp.reshape(overall_exp.shape[0], -1)))                   
            elif features == 1:
                feat_arr_fin = np.array([]).reshape(0, 64)
                for feature in range(len(comb_dicts[event][features])):
                    feat_arr = comb_dicts[event][features][feature].reshape(1, -1)
                    feat_arr_fin = np.vstack((feat_arr_fin, feat_arr))
                fin_arr = np.hstack((fin_arr, feat_arr_fin))
            else:
                feat_arr_fin = np.array([]).reshape(0, 1)
                for feature in range(len(comb_dicts[event][features])):
                    feat_arr = np.array(int(comb_dicts[event][features][feature]))
                    feat_arr_fin = np.vstack((feat_arr_fin, feat_arr))
                fin_arr = np.hstack((fin_arr, feat_arr_fin))
        fin_arr_all = np.vstack((fin_arr_all, fin_arr))
        
    return fin_arr_all
                