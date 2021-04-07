%% func l3 2-stage
function [A_l4_1stage,tot_l4_1stage] = Dic_Learn_l4(Y,A0)
[dic_size, sample_size] = size(Y);
A_l4_2stage = A0;
last_A_l4_2stage = zeros(size(A0));
e_l4_2stage = 1;
iter_l4_2stage = 1;
% f_val_l4_2stage = [];
% conv_l4_2stage = [];
tic
while (e_l4_2stage > 1e-4 && iter_l4_2stage < 50)
   
    AY_l4_2stage =  A_l4_2stage'*Y;
    phase_l4_2stage = 1./(abs(AY_l4_2stage)).*(AY_l4_2stage);
    phase_l4_2stage (find(isnan(phase_l4_2stage)))=0;
    dA_l4_2stage = (4 * abs(AY_l4_2stage).^3 .*    phase_l4_2stage * Y')';%AY.^(lp-1)*Y'/m;
   
   
    [U_l4_2stage,~,V_l4_2stage] = svd(-dA_l4_2stage,'econ');%-12*m*p^2*A%+2*p*noi_var^2+noi_var^4
    A_l4_2stage = U_l4_2stage*V_l4_2stage';
    
%     f_val_l4_2stage = [f_val_l4_2stage sum(sum((abs(A_l4_2stage' * Y))))/sample_size];
    e_l4_2stage = norm(A_l4_2stage-last_A_l4_2stage,'fro')/dic_size;
    
%     conv_l4_2stage = [conv_l4_2stage 1-sum(sum(abs(A_l4_2stage'*D0).^4))/dic_size];
    last_A_l4_2stage = A_l4_2stage;
    iter_l4_2stage =  iter_l4_2stage + 1;
    
end
A_l4_1stage=last_A_l4_2stage;
tot_l4_1stage = toc;