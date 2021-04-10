function [Error,Time]= TSP_rev1_syndata (N,K,M,Theta,SNR)
% SNR = 0;%norm(Y0,'fro')/sqrt(n*m)*sqrt(0);
% Theta_corr = 0.2;
Test_num = 50;
Error_l3_2stage = 0;
Error_l3_1stage = 0;
Error_l4_MSP = 0;
Error_l1_RGD = 0;

Error_logcosh_RTR = 0;
Error_KSVD = 0;
Error_SPAMS = 0;
Error_Trans = 0;

Time_l3_2stage = 0;
Time_l3_1stage = 0;
Time_l4_MSP = 0;
Time_l1_RGD = 0;

Time_logcosh_RTR = 0;
Time_KSVD = 0;
Time_SPAMS = 0;
Time_Trans = 0;
for trial_ind = 1:Test_num
    %% Gen dict
    %rng(2019)
    Dini = randn(N);
    [Q,~] = qr(Dini);
    D0 = Q(:,1:K);
    %% Gen sparse code
    B = rand(K,M);
    B(B>Theta) = 0;
    B(B>0) = 1;
    G = randn(K,M);
    X0 = B.*G;
    %% Gen Measure
    Y0 = D0*X0;
    %% add noise
    raw_power = norm(Y0,'fro')^2/numel(Y0);
    noi_var = raw_power/(10^(SNR/10));
    Y = Y0 + (eps+noi_var)*randn(N,M);
    
    %% DL initial
    [dic_size, sample_size] = size(Y);
    [A0,~] = qr(randn(dic_size));
    A0 = A0(:,1:dic_size);
    %% L3 DL
      fprintf(' DL proposed start!\n');
    [A_l3_2stage_2,tot_l3_2stage] = Dic_Learn_l3_2stage (Y,A0);
    A_l3_2stage_2=De_permutation(A_l3_2stage_2,D0);
    Error_l3_2stage_inst  = norm( A_l3_2stage_2  - D0,'fro')^2/norm(D0,'fro')^2;
  
    %% l3 GPM
    fprintf(' DL L3 GPM start!\n');
    [A_l3_GPM,tot_l3_1stage] = Dic_Learn_l3_1stage (Y,A0);
    A_l3_GPM=De_permutation(A_l3_GPM,D0);
    %error
    Error_l3_1stage_inst  = norm( A_l3_GPM  - D0,'fro')^2/norm(D0,'fro')^2;
    
    %% L4
    fprintf(' DL L4 MSP start!\n');
    [A_l4_MSP,tot_l4_MSP] = Dic_Learn_l4(Y,A0);
    A_l4_MSP=De_permutation(A_l4_MSP,D0);
    %error
    Error_l4_MSP_inst = norm( A_l4_MSP  - D0,'fro')^2/norm(D0,'fro')^2;
    
    %% L1 DL
       fprintf(' DL L1 RGD start!\n');
    [A_l1_RGD,tot_l1_RGD] = Dic_Learn_l1_sub (Y,D0); %D0 to stop the alg
    A_l1_RGD=De_permutation(A_l1_RGD,D0);
    %error
    Error_l1_RGD_inst = norm( A_l1_RGD  - D0,'fro')^2/norm(D0,'fro')^2;
 
    %% logcosh RGD

    
    
    %% logcosh RTR
    fprintf(' DL Logcosh RTR start!\n');
    [A_logcosh_rtr,tot_logcosh_rtr] = Dic_Learn_logcosh_RTR (Y,D0);%D0 to stop the alg
    A_logcosh_rtr=De_permutation(A_logcosh_rtr,D0); 
   
    Error_logcosh_rtr_inst = norm( A_logcosh_rtr  - D0,'fro')^2/norm(D0,'fro')^2;
    
    %% KSVD
fprintf(' DL KSVD start!\n');
    [A_KSVD,tot_KSVD] = Dic_Learn_KSVD (Y,A0,Theta);
    A_KSVD=De_permutation(A_KSVD,D0);
    %error
    Error_KSVD_inst= norm( A_KSVD  - D0,'fro')^2/norm(D0,'fro')^2;
    
    %% Spams
    fprintf(' DL SPAMS start!\n');
    [A_Spams,tot_Spams] = Dic_Learn_Spams (Y,A0);
    A_Spams=De_permutation(A_Spams,D0);
    Error_Spams_inst= norm( A_Spams  - D0,'fro')^2/norm(D0,'fro')^2;
    
%% Trans
fprintf(' DL TransL start!\n');
    [A_Dic_TransL,tot_Dic_TransL,~] = Dic_Learn_TransL (Y,A0,Theta);
    A_TransL=De_permutation(A_Dic_TransL,D0);
    %error
    Error_TransL_inst= norm( A_TransL - D0,'fro')^2/norm(D0,'fro')^2;
    
    %% Display
    fprintf('trail %d\n',  trial_ind );
    fprintf('SNR: %d ,  N: %d, K: %d, M: %d, theta: %d\n', SNR, N,K,M,Theta);
    fprintf('Proposed  error : %g, time %e\n', Error_l3_2stage_inst,tot_l3_2stage);
    fprintf('l3_GPM  error : %g, time %e\n', Error_l3_1stage_inst  , tot_l3_1stage);
    fprintf('l4_MSP  error : %g, time %e\n', Error_l4_MSP_inst ,tot_l4_MSP);
    fprintf('l1_RGD  error : %g, time %e\n', Error_l1_RGD_inst ,tot_l1_RGD);
    fprintf('LCH_RTR  error : %g, time %e\n', Error_logcosh_rtr_inst ,tot_logcosh_rtr);
    fprintf(' KSVD  error : %g, time %e\n', Error_KSVD_inst , tot_KSVD );
    fprintf(' SPAMS  error : %g, time %e\n', Error_Spams_inst , tot_Spams );
    fprintf(' Trans  error : %g, time %e\n',  Error_TransL_inst , tot_Dic_TransL );
    
    %% cal metric
   
    Error_l3_2stage = Error_l3_2stage  + Error_l3_2stage_inst/Test_num;
    Error_l3_1stage = Error_l3_1stage  +   Error_l3_1stage_inst/Test_num;
    Error_l4_MSP = Error_l4_MSP + Error_l4_MSP_inst/Test_num;
    Error_l1_RGD = Error_l1_RGD + Error_l1_RGD_inst/Test_num;
    Error_logcosh_RTR = Error_logcosh_RTR+Error_logcosh_rtr_inst/Test_num;
    Error_SPAMS= Error_SPAMS + Error_Spams_inst/Test_num;
    Error_KSVD= Error_KSVD + Error_KSVD_inst/Test_num;
    Error_Trans = Error_Trans+Error_TransL_inst/Test_num;
    Time_l3_2stage = Time_l3_2stage + tot_l3_2stage/Test_num;
    Time_l3_1stage = Time_l3_1stage + tot_l3_1stage/Test_num;
    Time_l4_MSP = Time_l4_MSP + tot_l4_MSP/Test_num;
    Time_l1_RGD = Time_l1_RGD+ tot_l1_RGD/Test_num;
    Time_logcosh_RTR = Time_logcosh_RTR+tot_logcosh_rtr/Test_num;
    Time_KSVD = Time_KSVD + tot_KSVD/Test_num;
    Time_SPAMS = Time_SPAMS + tot_Spams/Test_num;
    Time_Trans = Time_Trans+tot_Dic_TransL/Test_num;
end
Error = [Error_l3_2stage, Error_l3_1stage, Error_l4_MSP,Error_l1_RGD,...
    Error_logcosh_RTR,Error_KSVD,Error_SPAMS,Error_Trans];
Time = [Time_l3_2stage, Time_l3_1stage ,  Time_l4_MSP, Time_l1_RGD,...
    Time_logcosh_RTR,Time_KSVD,Time_SPAMS,Time_Trans];
