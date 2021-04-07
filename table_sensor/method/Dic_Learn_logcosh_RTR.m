function[A_logcosh_rtr,tot_logcosh_rtr] = Dic_Learn_logcosh_RTR (Y,D0)
%% func logcosh

    [n,~]= size(D0);
    [~,m]= size(Y);
    mu =1/100;
    tic

    Q = [];

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
        q = Proj_TR_Sphere( Y, mu, q_i,P,D0);
        r = P*q;
        rp= rounding_sphere( r,Y);
       
        
%         rp = l1_rounding(r,Y);
        Q = [Q rp];
%          sg=  sum(abs(rp'*D0).^4)

end
    rp = null(Q')/norm(null(Q'));
    Q = [Q rp];
A_logcosh_rtr = Q;

tot_logcosh_rtr = toc;