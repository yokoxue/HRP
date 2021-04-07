%% func l3 2-stage
function [A_l3_2stage_2,tot_l3_2stage] = Dic_Learn_l3_2stage (Y,A0)
[dic_size, sample_size] = size(Y);
A_l3_2stage = A0;
last_A_l3_2stage = zeros(size(A0));
e_l3_2stage = 1;
iter_l3_2stage = 1;
% f_val_l3_2stage = [];
% conv_l3_2stage = [];
tic
while (e_l3_2stage > 1e-4 && iter_l3_2stage < 50)
% tic
    
    AY_l3_2stage =  A_l3_2stage'*Y;
%     phase_l3_2stage = 1./(abs(AY_l3_2stage)).*(AY_l3_2stage);
%      phase_l3_2stage (find(isnan(phase_l3_2stage)))=0;
    dA_l3_2stage = (3 * abs(AY_l3_2stage) .*    AY_l3_2stage  * Y')';%AY.^(lp-1)*Y'/m;
   
    [U_l3_2stage,~,V_l3_2stage] = svd(-dA_l3_2stage,'econ');%-12*m*p^2*A%+2*p*noi_var^2+noi_var^4
    A_l3_2stage = U_l3_2stage*V_l3_2stage';
    
%     f_val_l3_2stage = [f_val_l3_2stage sum(sum((abs(A_l3_2stage' * Y).^4)))/sample_size];
    e_l3_2stage = norm(A_l3_2stage-last_A_l3_2stage,'fro')/dic_size;
    
%     conv_l3_2stage = [conv_l3_2stage 1-sum(sum(abs(A_l3_2stage'*D0).^4))/dic_size];
    last_A_l3_2stage = A_l3_2stage;
    iter_l3_2stage =  iter_l3_2stage + 1;
%    tl3 = toc
end

r = last_A_l3_2stage;

[n,p] = size(Y);
% Obj = @(h) 1/(p)  *sum(sum(abs( h'*Y_p)));
q_0 = r; %initialization
q = q_0;
mu = 0.1; rho = 0.8;  mu_threshold = 1e-10;
conv = [];
iter = 1;
while (iter < 50&&mu> mu_threshold)
    grad = 1/(p) *  Y*sign(Y'*q);
    mu = mu*rho;
  
    q_new = r + 0.5*((q - mu*grad)-r*(q - mu*grad)'*r);
  
    q = q_new; %update
    iter = iter +1;
   
end
[U_l3_2stage,~,V_l3_2stage] = svd(q,'econ');%-12*m*p^2*A%+2*p*noi_var^2+noi_var^4
A_l3_2stage_2 = U_l3_2stage*V_l3_2stage';
% [A_l3_2stage_2, conv_l3_2stage_2]=rounding_new_st(Y,last_A_l3_2stage,zeros(size(A_l3_2stage)));
tot_l3_2stage = toc;
