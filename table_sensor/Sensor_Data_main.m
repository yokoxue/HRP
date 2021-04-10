basePath = [fileparts(mfilename('fullpath')) filesep];
addpath(genpath(basePath ));
%import excel
Pm_25 = load('poland_pm.mat');
Pm_25 = (Pm_25.Pm_25);
%ratio = 1;
%numer_sample = size(Pm_25,1)*ratio;
pm25_Raw = Pm_25(1:744,:);
pm25_miss_ind = find(isnan(pm25_Raw));
compset = setdiff([1:1:numel(pm25_Raw)]',pm25_miss_ind);

% deal missing data
pm25_Raw_comp = zeros(size(pm25_Raw ));
for ti = 1:1:size(pm25_Raw,1)
    temp = pm25_Raw(ti,:);
    temp( find (isnan(temp))) = randn(size( find (isnan(temp))))+mean(find (~isnan(temp)));
    %      temp( find (isnan(temp))) = mean(find (~isnan(temp)));
    % temp( find (isnan(temp))) = 0;
    pm25_Raw_comp (ti,:) =  temp;
end
%normalize
pm25_data_Max = max(max(pm25_Raw_comp));
pm25_Normalized = pm25_Raw_comp./pm25_data_Max;
mean25 = mean(pm25_Normalized);
numer_sample = size(pm25_Raw,1);
Y = [pm25_Normalized(:,:)'];


%% DL
[dic_size, sample_size] = size(Y);
[A0,~] = qr(randn(dic_size));
A0 = A0(:,1:dic_size);
%% prop DL
[A_l3,tot_l3] = Dic_Learn_l3_2stage (Y,A0);
%% L3
[A_l3_s1,tot_l3_s1] = Dic_Learn_l3_1stage (Y,A0);
%% L4
[A_l4,tot_l4] = Dic_Learn_l4 (Y,A0);
%% L1 DL
[A_l1,tot_l1_FIX] = l1_O_RGD  (Y,A_l3);
%% RTR
[A_logcosh_rtr,tot_logcosh_rtr] = Dic_Learn_logcosh_RTR (Y,A_l3);
%% Spams
[A_Spams,tot_Spams] = Dic_Learn_Spams (Y,A0);



%%
%% Sparse coding
Sparse_code_l3 = A_l3'*Y ;
Sparse_code_l3_s1 = A_l3_s1'*Y ;
Sparse_code_l4 = A_l4'*Y ;
Sparse_code_l1 = A_l1'*Y ;
Sparse_code_log = A_logcosh_rtr'*Y ;
Sparse_code_SPAMS = inv(A_Spams)*Y ;

%% compress
TIME_l3 = tot_l3*ones(1,size(A_l3,1));
TIME_l3_s1 = tot_l3_s1*ones(1,size(A_l3,1));
TIME_l4 = tot_l4*ones(1,size(A_l3,1));
TIME_l1 = tot_l1_FIX*ones(1,size(A_l3,1));
TIME_log = ones(1,size(A_l3,1));
TIME_SPAMS= ones(1,size(A_l3,1));



RMSE_l3 = zeros(1,size(A_l3,1));
RMSE_l3_s1 = zeros(1,size(A_l3,1));
RMSE_l4 = zeros(1,size(A_l3,1));
RMSE_l1 = zeros(1,size(A_l3,1));
RMSE_log = zeros(1,size(A_l3,1));
RMSE_SPAMS= zeros(1,size(A_l3,1));
RMSE_KSVD= zeros(1,size(A_l3,1));
RMSE_TransL= zeros(1,size(A_l3,1));
TIME_KSVD = zeros(1,size(A_l3,1));
TIME_TransL = zeros(1,size(A_l3,1));
for num_b = [5,7,11,18,28]
    
    disp('Compressed');disp(num_b);
    pm25_recon_sample_l3  = zeros(size(Y ));
    pm25_recon_sample_l3_s1  = zeros(size(Y ));
    pm25_recon_sample_l4  = zeros(size(Y ));
    pm25_recon_sample_l1  = zeros(size(Y ));
    pm25_recon_sample_log  = zeros(size(Y ));
    pm25_recon_sample_SPAMS  = zeros(size(Y ));
    
    for ci = 1:1:size(Sparse_code_l3,2)
        pm25_recon_sample_l3(:,ci) = Recon_col(Sparse_code_l3(:,ci),A_l3,num_b);
        pm25_recon_sample_l3_s1(:,ci) = Recon_col(Sparse_code_l3_s1(:,ci),A_l3_s1,num_b);
        pm25_recon_sample_l4(:,ci) = Recon_col(Sparse_code_l4(:,ci),A_l4,num_b);
        pm25_recon_sample_l1(:,ci) = Recon_col(Sparse_code_l1(:,ci),A_l1,num_b);
        pm25_recon_sample_log(:,ci) = Recon_col(Sparse_code_log(:,ci),A_logcosh_rtr,num_b);
    
     pm25_recon_sample_SPAMS(:,ci)=Recon_col(Sparse_code_SPAMS(:,ci), A_Spams,num_b);
    end
    %% KSVD
    [A_KSVD,X_KSVD,tot_KSVD] = Dic_Learn_KSVD (Y,A0,num_b/dic_size);
    %% Trans
    [A_TransL,tot_TransL,X_TransL]  =Dic_Learn_TransL (Y,A0,num_b/dic_size);
    %% error
    pm25_recon_sample_l3 = pm25_recon_sample_l3'.*pm25_data_Max;
    pm25_recon_sample_l3_s1 = pm25_recon_sample_l3_s1'.*pm25_data_Max;
    pm25_recon_sample_l4 = pm25_recon_sample_l4'.*pm25_data_Max;
    pm25_recon_sample_log = pm25_recon_sample_log'.*pm25_data_Max;
    pm25_recon_sample_l1 = pm25_recon_sample_l1'.*pm25_data_Max;
    pm25_recon_sample_KSVD = (A_KSVD*X_KSVD)'.*pm25_data_Max;
   pm25_recon_sample_SPAMS = pm25_recon_sample_SPAMS'.*pm25_data_Max;
    pm25_recon_sample_TransL = (A_TransL*X_TransL)'.*pm25_data_Max;
    
    RMSE_l3(1,num_b ) = sqrt(norm(pm25_recon_sample_l3(compset)-pm25_Raw_comp(compset),'fro')^2/(norm(pm25_Raw_comp(compset),'fro'))^2);
    RMSE_l3_s1(1,num_b ) = sqrt(norm(pm25_recon_sample_l3_s1(compset)-pm25_Raw_comp(compset),'fro')^2/(norm(pm25_Raw_comp(compset),'fro'))^2);
    RMSE_l4(1,num_b ) = sqrt(norm(pm25_recon_sample_l4(compset)-pm25_Raw_comp(compset),'fro')^2/(norm(pm25_Raw_comp(compset),'fro'))^2);
    RMSE_log(1,num_b ) = sqrt(norm(pm25_recon_sample_log(compset)-pm25_Raw_comp(compset),'fro')^2/(norm(pm25_Raw_comp(compset),'fro'))^2);
    RMSE_l1(1,num_b ) = sqrt(norm(pm25_recon_sample_l1(compset)-pm25_Raw_comp(compset),'fro')^2/(norm(pm25_Raw_comp(compset),'fro'))^2);
    RMSE_KSVD(1,num_b ) = sqrt(norm(pm25_recon_sample_KSVD(compset)-pm25_Raw_comp(compset),'fro')^2/(norm(pm25_Raw_comp(compset),'fro'))^2);
    RMSE_SPAMS(1,num_b ) = sqrt(norm(pm25_recon_sample_SPAMS(compset)-pm25_Raw_comp(compset),'fro')^2/(norm(pm25_Raw_comp(compset),'fro'))^2);
    RMSE_TransL(1,num_b ) = sqrt(norm(pm25_recon_sample_TransL(compset)-pm25_Raw_comp(compset),'fro')^2/(norm(pm25_Raw_comp(compset),'fro'))^2);
   
    TIME_KSVD(1,num_b) = tot_KSVD;
    TIME_TransL(1,num_b) = tot_TransL;
end
