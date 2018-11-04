function y3=alg_wgrmf(Y,Sd,St,cv_setting,nr_fold,left_out,use_WKNKN,K,eta,use_W_matrix)
%alg_wgrmf predicts DTIs based on the WGRMF algorithm described in the following paper: 
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

    % received parameters are passed to 'alg_grmf'
    y3=alg_grmf(Y,Sd,St,cv_setting,nr_fold,left_out,use_WKNKN,K,eta,use_W_matrix);

    % behavior inside 'alg_grmf.m', 'alg_grmf_predict.m' & 'alg_grmf_parest.m'
    % changes based on the value of use_W_matrix, which in this case, is 1.

    % also note the value of use_W_matrix when the classifier is set as
    % 'wgrmf' in Run_Cross_validation.m

end