function [A_TransL,tot_TransL,X_TransL]  =Dic_Learn_TransL (Y,A0,Theta);
[N,~]=size(A0);
W =eye(N);
s= floor(Theta*N);
lambda = 0.5;
mu = 0.25;
eta = 10^-4;Tr_Iter=60;
Iter = 1000;
conv = [];
tic
for i = 1:1:Iter
    X=Trans_Sparsecode(W,Y,s);
    W=Trans_Dic(W,Y,X,lambda,mu,eta,Tr_Iter);
  %  obj  = Trans_Learn_Obj(Y,W,X,lambda,mu);
   % conv =[conv obj];
end
X_TransL = X;
%Xt = W*Y;
Wt = normc(W.');
A_TransL = Wt;
tot_TransL = toc;
%Dic_out=Y*Xt' *inv(Xt*Xt');
