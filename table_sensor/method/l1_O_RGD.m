%% L1 RGD Complete
function [A_l1_FIX,tot_l1_FIX] = l1_O_RGD (Y,D0)
[n,~]= size(D0);
[~,m]= size(Y);
 Q = [];
tic
    for i = 1:n-1
        if(i>1)
            P = null(Q');
            %         Y = P*Y;
        else
            P = eye(n);
        end
        q_i = randn(n-i+1,1);
        % %     q_i = P*q_i;
        q_i = q_i / norm(q_i);
        %     %     q = TR_Sphere(Y,mu,q_i,P,D0);
        q = l1_RGDcxxor(Y,D0,P,q_i,m);
        r = P*q;
        rp =r;
%         rp = l1_rounding(r,Y);
        Q = [Q rp];
%          sg=  sum(abs(rp'*D0).^4)

    end
    rp = null(Q')/norm(null(Q'));
    Q = [Q rp];
    A_l1_FIX = Q;
    tot_l1_FIX = toc;
%     sg=  sum((rp'*D0).^4)
%     H0 = [H0;rp'*Yr];
% 