function y3=alg_cmf(Y,Sd,St,cv_setting,nr_fold,left_out,use_WKNKN,K,eta,use_W_matrix)
%alg_cmf predicts DTIs based on the algorithm described in the following paper: 
% Xiaodong Zheng, Hao Ding, Hiroshi Mamitsuka and Shanfeng Zhu
% (2013) Collaborative Matrix Factorization with Multiple Similarities for Predicting Drug-Target Interactions
%
% INPUT:
%  Y:           interaction matrix
%  Sd:          pairwise drug similarities matrix
%  St:          pairwise target similarities matrix
%  cv_setting:  cross validation setting ('cv_d', 'cv_t' or 'cv_p')
%  nr_fold:     number of folds in cross validation experiment
%  left_out:    if cv_setting=='cv_d' --> left_out is 'drug' indices that are left out
%               if cv_setting=='cv_t' --> left_out is 'target' indices that are left out
%               if cv_setting=='cv_p' --> left_out is 'drug-target pair' indices that are left out
%
% OUTPUT:
%  y3:  prediction matrix

    % get best parameters
    [k,lambda_l,lambda_d,lambda_t,num_iter] = alg_cmf_parest(Y,Sd,St,cv_setting,nr_fold,left_out,use_WKNKN,K,eta,use_W_matrix);
    fprintf('k%g\t\t%g\t%g\t%g\t\t',k,lambda_l,lambda_d,lambda_t);

    % preprocessing Y
    if use_WKNKN
        Y = preprocess_WKNKN(Y,Sd,St,K,eta);
    end

    % initialize A & B
    [A,B] = initializer(Y,k);

    % predict
    test_ind = get_test_indices(Y,cv_setting,left_out);
    W = ones(size(Y));
    W(test_ind) = 0;
    [A,B] = alg_cmf_predict(Y,A,B,Sd,St,lambda_l,lambda_d,lambda_t,num_iter,W);

    % compute prediction matrix
    y3 = A*B';

    %--------------------------------------------------------------------

end