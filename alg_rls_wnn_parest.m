function eta=alg_rls_wnn_parest(Y,Sd,St,cv_setting,nr_fold,left_out)
%alg_rls_wnn_parest is a helper function of RLS-WNN that estimates the
%parameter, eta.

    %--------------------------------------------------------------------

    % ranges of parameter values to be tested to identify best combination
    eta_range = 0.05:0.05:1;

    %--------------------------------------------------------------------

    folds = get_folds(Y,cv_setting,nr_fold,left_out);   % folds of the CV done on the training set

    %--------------------------------------------------------------------

    y2s = cell(1,nr_fold);
    for i=1:nr_fold
        y2 = Y;
        y2(folds{i}) = 0;  % folds{i} is the validation set
        y2s{i} = y2;
    end

    %--------------------------------------------------------------------

    best_AUPR = -Inf;
    for t_indx=1:length(eta_range)
        eta = eta_range(t_indx);
        
        AUPRs = zeros(nr_fold,1);
        for z=1:nr_fold

            % PREDICT
            y2 = y2s{z};

            % WNN -----------------------------------------------------
            if strcmp(cv_setting,'cv_d')
                empty_rows = find(any(Y,2) == 0);   % empty rows
                w = eta .^ (0:length(Sd)-1);
                w(w < 10^-4) = [];
                k = length(w);
                for r=1:length(empty_rows)
                    i = empty_rows(r);
                    drug_sim = Sd(i,:);
                    drug_sim(i) = 0;    % set self-similarity to ZERO
                    [~,indx]=sort(drug_sim,'descend');
                    indx = indx(1:k);
                    y2(i,:) = w * Y(indx,:);
                end
            elseif strcmp(cv_setting,'cv_t')
                empty_cols = find(any(Y) == 0);   % empty columns
                w = eta .^ (0:length(St)-1);
                w(w < 10^-4) = [];
                k = length(w);
                for c=1:length(empty_cols)
                    j = empty_cols(c);
                    target_sim = St(j,:);
                    target_sim(j) = 0;    % set self-similarity to ZERO
                    [~,indx]=sort(target_sim,'descend');
                    indx = indx(1:k);
                    y2(:,j) = Y(:,indx) * w';
                end
            end
            % ---------------------------------------------------------

            % GIP -----------------------------------------------------
            alpha = 0.5;
            Sd = alpha*Sd + (1-alpha)*getGipKernel(y2);
            St = alpha*St + (1-alpha)*getGipKernel(y2');
            % ---------------------------------------------------------

            sigma = 1;

            [va,la] = eig(Sd);
            [vb,lb] = eig(St);

            l = kron(diag(lb)',diag(la));
            l = l ./ (l + sigma);

            m1 = va' * y2 * vb;
            m2 = m1 .* l;
            y3 = va * m2 * vb';

            [~,AUPRs(z)] = returnEvaluationMetrics(Y(folds{z}),y3(folds{z}));
        end
        aupr_res = mean(AUPRs);
        if best_AUPR < aupr_res
            best_AUPR = aupr_res;
            best_eta = eta;
        end
    end

    %--------------------------------------------------------------------

    % return best parameters
    eta = best_eta;

    %--------------------------------------------------------------------
    
end