def ml(feat_and_lab_txt, tmp_dir, num_rna_seq_samples):
    
    X_train = np.array([]).reshape(0, 81).astype('float32')
    with gzip.open(tmp_dir+feat_and_lab_txt, 'rt') as f:
        for line in f:
            key, dinuc, rna_exp, ssv, wsv, ca, lab = line.strip('/n').split('\t')
            dinuc_list, rna_exp_list, ssv_list, wsv_list, ca_list, lab_list = dinuc.split(), rna_exp.split(), ssv.split(), wsv.split(), ca.split(), lab.split() 
            dinuc_arr, pre_rna_exp_arr, rna_exp_arrs, rna_exp_arr, ssv_arr, wsv_arr, ca_arr, lab_arr   = [], [], [], [], [], [], [], []
            trans_size = int((len(rna_exp_list)-num_rna_seq_samples)/num_rna_seq_samples)
            for i in range(len(dinuc_list)):
                dinuc_arr.append(np.array([float(x) for x in dinuc_list[i].split(',')])) 
            for i in range(len(rna_exp_list)):
                pre_rna_exp_arr.append(np.array([float(x) for x in rna_exp_list[i].split(',')])) 
            for samps in range(1, num_rna_seq_samples+1):
                rna_exp_arrs.append(pre_rna_exp_arr[int(samps*trans_size-trans_size+samps-1):int(samps*trans_size+samps)])
            elem_comb_rna_exp_arrs = list(zip(*rna_exp_arrs))
            rna_exp_arr = [np.mean(x, axis = 0) for x in elem_comb_rna_exp_arrs]
            for i in range(len(ssv_list)):
                ssv_arr.append(np.array([float(x) for x in ssv_list[i].split(',')])) 
            for i in range(len(wsv_list)):
                wsv_arr.append(np.array([float(x) for x in wsv_list[i].split(',')])) 
            for i in range(len(ca_list)):
                ca_arr.append(np.array([float(x) for x in ca_list[i].split(',')])) 
            for i in range(len(lab_list)):
                lab_arr.append(np.array([float(x) for x in lab_list[i].split(',')]))
            for i in list(zip(dinuc_arr, rna_exp_arr[:-1], np.repeat(rna_exp_arr[-1], len(dinuc_list)),  ssv_arr, wsv_arr, ca_arr, lab_arr)):
                ind_feat = np.hstack(i).reshape(1, -1).astype('float32')
                X_train = np.vstack((X_train, ind_feat))       
        
    X_train, X_test, y_train, y_test = train_test_split(feat_and_lab_arr[:, 0:80], feat_and_lab_arr[:, -1].astype('int32'), random_state = 0)
    
    scaler = MinMaxScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)
    
    smote = SMOTE(random_state = 0)
    X_train, y_train = smote().fit_sample(X_train, y_train)
    X_train, empty_1, y_train, empty_2 = train_test_split(X_train, y_train, random_state = 0, test_size = 0)
    
    c_range = [0.1, 0.5, 1, 10, 100]
    max_iter_range = [10000, 100000]
    solver_options = ['liblinear', 'saga']
    penalty_options = ['l1', 'l2']
    param_grid = dict(C = c_range, max_iter = max_iter_range, solver = solver_options, penalty = penalty_options)
    
    f_two_scorer = make_scorer(fbeta_score, beta = 2)
    f_half_scorer = make_scorer(fbeta_score, beta = 0.5)
    scoring = ['f1', f_two_scorer, f_half_scorer, 'roc_auc', 'average_precision']
    
    lr = LogisticRegression(random_state = 0)
    
    fin_dict = {}
    for score in scoring:
        lr_grid = GridSearchCV(lr, param_grid = param_grid, scoring = score, cv = 5)
        lr_grid.fit(X_train, y_train)
        fin_dict['%s_optimized_smote_model' % (score)] = lr_grid
    fin_dict['upsampled_smote_data'] = ([X_train, X_test, y_train, y_test])
    fin_dict['scale'] = scaler
    
    return fin_dict
   
