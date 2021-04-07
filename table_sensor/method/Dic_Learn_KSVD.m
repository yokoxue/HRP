function[A_KSVD,X,tot_KSVD] = Dic_Learn_KSVD (Y,A0,Theta)

%% DL
N = size(Y,1);
params.data = Y;
params.Tdata = Theta*N ;
params.initdict = A0;
tic
[A_KSVD,X] = ksvd(params);
tot_KSVD = toc;

