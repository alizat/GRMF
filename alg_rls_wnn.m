function y3=alg_rls_wnn(Y,Sd,St,cv_setting,nr_fold,left_out,use_WKNKN,K,eta,~)
%alg_rls_wnn predicts DTIs based on the algorithm described in the following paper: 
% Twan van Laarhoven and Elena Marchiori
% (2013) Predicting Drug-Target Interactions for New Drug Compounds Using a Weighted Nearest Neighbor Profile
%
% Modified from code of:
%  Twan van Laarhoven, Sander B. Nabuurs, Elena Marchiori,
%  (2011) Gaussian interaction profile kernels for predicting drug–target interaction
%  http://cs.ru.nl/~tvanlaarhoven/drugtarget2013/
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

    %--------------------------------------------------------------------

    % parameters
    alpha = 0.5;
    %eta = alg_rls_wnn_parest(Y,Sd,St,cv_setting,nr_fold,left_out);

    %--------------------------------------------------------------------

    % preprocessing Y
    if use_WKNKN
        Y = preprocess_WKNKN(Y,Sd,St,K,eta);
    else
        Y = preprocess_WNN(Y,Sd,St,cv_setting,eta);
    end

    % GIP 
    Sd = alpha*Sd + (1-alpha)*getGipKernel(Y);
    St = alpha*St + (1-alpha)*getGipKernel(Y');


    % RLS_KRON
	sigma = 1;
	
	[va,la] = eig(Sd);
	[vb,lb] = eig(St);
	
	l = kron(diag(lb)',diag(la));
	l = l ./ (l + sigma);
	
	m1 = va' * Y * vb;
	m2 = m1 .* l;
	y3 = va * m2 * vb';

    %fprintf('eta=%g\t\t',eta);

end