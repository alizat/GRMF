function [A,B]=alg_cmf_predict(Y,A,B,Sd,St,lambda_l,lambda_d,lambda_t,num_iter,W)
%alg_cmf_predict is a helper function of CMF that executes the alternating
%least squares method to obtain the solution. 
    
    K = size(A,2);
    lambda_d_Sd = lambda_d*Sd;
    lambda_t_St = lambda_t*St;
    lambda_l_eye_K = lambda_l*eye(K);

    if nargin < 10
        AtA = A'*A;
        BtB = B'*B;
        for z=1:num_iter
            A = (Y*B + lambda_d_Sd*A)  / (BtB + lambda_l_eye_K + lambda_d*(AtA));
            AtA = A'*A;
            B = (Y'*A + lambda_t_St*B) / (AtA + lambda_l_eye_K + lambda_t*(BtB));
            BtB = B'*B;
        end
        
    else
        H = W .* Y;
        for z=1:num_iter
            A_old = A;
            HB_plus_lambda_d_Sd_A_old = H*B + lambda_d_Sd*A_old;
            lambda_l_eye_k_plus_lambda_d_A_oldt_A_old = lambda_l_eye_K + lambda_d*(A_old'*A_old);
            for a=1:size(A,1)
                A(a,:) = HB_plus_lambda_d_Sd_A_old(a,:) / (B'*diag(W(a,:))*B + lambda_l_eye_k_plus_lambda_d_A_oldt_A_old);
            end
            B_old = B;
            HtA_plus_lambda_t_St_B_old = H'*A + lambda_t_St*B_old;
            lambda_l_eye_k_plus_lambda_t_B_oldt_B_old = lambda_l_eye_K + lambda_t*(B_old'*B_old);
            for b=1:size(B,1)
                B(b,:) = HtA_plus_lambda_t_St_B_old(b,:) / (A'*diag(W(:,b))*A + lambda_l_eye_k_plus_lambda_t_B_oldt_B_old);
            end

%             % for readability...
%             A_old = A;
%             lambda_d_A_oldt_A_old = lambda_d*(A_old'*A_old);
%             for a=1:size(A,1)
%                 A(a,:) = (H(a,:)*B + lambda_d_Sd(a,:)*A_old) / (B'*B + lambda_l_eye_k + lambda_d_A_oldt_A_old);
%             end
%             B_old = B;
%             lambda_t_B_oldt_B_old = lambda_t*(B_old'*B_old);
%             for b=1:size(B,1)
%                 B(b,:) = (H(:,b)'*A + lambda_t_St(b,:)*B_old) / (A'*A + lambda_l_eye_k + lambda_t_B_oldt_B_old);
%             end
        end
    end
    
end