function [A_De]=De_permutation(A_est,D0)
p_mat = A_est' *  D0;
sign_p = sign (p_mat);
[~,flagone] = max (abs(p_mat));
p_mat_0 = eye(size(p_mat));
p_mat_0 = p_mat_0 (:,flagone);
p_mat_1 = p_mat_0.*sign_p;
A_De= A_est * p_mat_1;