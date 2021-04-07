%% l1 RGD
function q1 = l1_RGDcxxor(Y,D0,U,q,m)
%             record_basis = zeros(n, 1);

% random initialization on sphere
eps = 10^-3;
max_iter =4000;
for i = 1:max_iter
%     Tq = null( (U*q)'); 
    eta = 1/sqrt(i);  % diminishing step size
    egrad = 1/m*U'*Y*sign(Y'*U*q);
    q = q - eta * (egrad - q*((q)'*egrad));    % Riemannian step
    q = q/norm(q);   % projection
    q_rot = D0' *U*q;    % should be 1-sparse if one dictionary element is found   
    % early stopping if recovered
    if abs(max(abs(q_rot))-1) <= eps
        [~,ind] = max(abs(q_rot));
        %record_basis(ind) = 1;  % recovery up to sign
        break
    end
    
end
q1=q;
% fprintf('Iterl1: %g\n',i);