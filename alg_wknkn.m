function y3=alg_wknkn(Y,Sd,St,~,~,~,~,K,eta,~)
%alg_wknkn predicts DTIs based on the WKNKN algorithm described in the following paper: 
% Ali Ezzat, Peilin Zhao, Min Wu, Xiao-Li Li and Chee-Keong Kwoh
% (2016) Drug-target interaction prediction with graph-regularized matrix factorization
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
    y3 = preprocess_WKNKN(Y,Sd,St,K,eta);
end