def ml(feat_and_lab_arr):
    
    X_train_raw, X_test_raw, y_train_raw, y_test = train_test_split(feat_and_lab_arr[:, 0:80], feat_and_lab_arr[:, -1].astype('int32'), random_state = 0)
    
    scaler = MinMaxScaler()
    X_train_raw_scaled = scaler.fit_transform(X_train_raw)
    X_test = scaler.transform(X_test_raw)
    
    smote = SMOTE(random_state = 0)
    X_train_unshuffled, y_train_unshuffled = smote().fit_sample(X_train_raw_scaled, y_train_raw)
    X_train, empty_1, y_train, empty_2 = train_test_split(X_train_unshuffled, y_train_unshuffled, random_state = 0, test_size = 0)
    
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
    
    return fin_dict
   
