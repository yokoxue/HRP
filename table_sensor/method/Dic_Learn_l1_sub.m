%% func l1 lixiao
function [A_l1_FIX,tot_l1_FIX] = Dic_Learn_l1_sub (Y,D0)
[n,~]= size(D0);
[~,m]= size(Y);
A_l1_FIX = zeros(size(D0));
record_basis = zeros(n, 1);
max_iter = 1500;
eps = 10^-2; 
tic
for l = 1:round(5*n*log(n))
    % random initialization on sphere
    q = randn(n,1);
    q = q./norm(q);   
    for i = 1:max_iter
        eta = 1/sqrt(i);  % diminishing step size
        egrad = 1/m*Y*sign(Y'*q);
        q = q - eta * (egrad - q*(q'*egrad));    % Riemannian step
        q = q/norm(q);   % projection   
        q_rot = D0' *q;    % should be 1-sparse if one dictionary element is found
        % early stopping if recovered
        t=abs(max(abs(q_rot))-1);
        if t <= eps
            [~,ind] = max(abs(q_rot));
            record_basis(ind) = 1;  % recovery up to sign
            A_l1_FIX(:,ind) = q;
            break
        end     
    end
    
    if nnz(~(record_basis)) == 0   % found all standard basis
        break
    end
    
end  % end for 5*n*log(n) runs
tot_l1_FIX = toc;