%% Main comp for TNNLS syn data
clear all; clc; warning('off');
basePath = [fileparts(mfilename('fullpath')) filesep];
addpath(genpath(basePath ));
 dic_size = [20];%for performance
 sample_size = flip(linspace(100,1000,10)); %for performance

theta = flip([0.1:0.1:0.7]); %for performance

% dic_size = [10:10:50];%for complex
 %sample_size = 1000;%for complex
% 
theta = 0.2;%for complex
% N = dic_size;
% K = N;
SNR = inf;
Theta_corr = 0;
Num_of_method = 8;
Error = zeros(length(sample_size),length(dic_size ),Num_of_method );
Time = zeros(length(sample_size),length(dic_size ),Num_of_method );
for si = 1:1:length(sample_size)
    for ti = 1:1:length(dic_size )
        Theta = theta;
        M = sample_size(si);
        N = dic_size(ti);
        K=N;
        [Error(si,ti,:),Time(si,ti,:)]= Syndata_test(N,K,M,Theta,SNR);
    end
end

% save('Error_samp_vs_dic0918.mat','Error');
% save('Comp_Err210330_dic101050sample5the02.mat','Error');
 %save('Comp_Time210330_dic101050samplepaperthe02.mat','Time');
