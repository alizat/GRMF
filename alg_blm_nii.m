function y3=alg_blm_nii(Y,Sd,St,cv_setting,~,left_out,use_WKNKN,K,eta,~)
%alg_blm_nii predicts DTIs based on the algorithm described in the following paper: 
% Jian-Ping Mei, Chee-Keong Kwoh, Peng Yang, Xiao-Li Li and Jie Zheng
% (2013) Drug-target interaction prediction by learning from local information and neighbors
%
% INPUT:
%  Y:           interaction matrix
%  Sd:          pairwise drug similarities matrix
%  St:          pairwise target similarities matrix
%  cv_setting:  cross validation setting ('cv_d', 'cv_t' or 'cv_p')
%  left_out:    if cv_setting=='cv_d' --> left_out is 'drug' indices that are left out
%               if cv_setting=='cv_t' --> left_out is 'target' indices that are left out
%               if cv_setting=='cv_p' --> left_out is 'drug-target pair' indices that are left out
%
% OUTPUT:
%  y3:  prediction matrix

    %--------------------------------------------------------------------

    % parameters
    alpha = 0.5;
    delta = 1;

    %--------------------------------------------------------------------

    % get test indices
    test_ind = get_test_indices(Y,cv_setting,left_out);

    % preprocessing Y
    if use_WKNKN
        Y = preprocess_WKNKN(Y,Sd,St,K,eta);
    end


    % BLM-NII
    y3d = zeros(size(Y));
    y3t = zeros(size(Y));
    for k=1:length(test_ind)
        int = test_ind(k);
        [d,t] = ind2sub(size(Y),int);


        % drug side prediction ----------------------------------
        y_temp = Y;
        % NII ---------------------------------
        if sum(Y(d,:))==0
            y_temp(d,:) = Sd(d,:) * Y;
            maxx = max(y_temp(d,:));
            minn = min(y_temp(d,:));
            y_temp(d,:) = (y_temp(d,:) - minn) / (maxx - minn); % normalization
        end
        % end NII -----------------------------
        train = 1:length(y_temp(d,:));   train(t) = [];
        test = t;
        Kt = alpha*St + (1-alpha)*getGipKernel(y_temp');
        K1 = Kt(train,train);
        K2 = Kt(test,train);
        n = length(train);
        model = (K1+delta*eye(n)) \ y_temp(d,train)';
        y3d(d,t) = K2*model;
        % end drug side prediction ------------------------------


        % target side prediction --------------------------------
        y_temp = Y;
        % NII ---------------------------------
        if sum(Y(:,t))==0
            y_temp(:,t) = Y * St(t,:)';
            maxx = max(y_temp(:,t));
            minn = min(y_temp(:,t));
            y_temp(:,t) = (y_temp(:,t) - minn) / (maxx - minn); % normalization
        end
        % end NII -----------------------------
        train = 1:length(y_temp(:,t));   train(d) = [];
        test = d;
        Kd = alpha*Sd + (1-alpha)*getGipKernel(y_temp);
        K1 = Kd(train,train);
        K2 = Kd(test,train);
        n = length(train);
        model = (K1+delta*eye(n)) \ y_temp(train,t);
        y3t(d,t) = K2*model;
        % end target side prediction ----------------------------
    end

    y3 = (y3d + y3t) / 2;   % g=avg
    %y3 = max(y3d,y3t);      % g=max

end