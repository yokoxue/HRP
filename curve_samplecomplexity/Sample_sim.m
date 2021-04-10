
% SNR = 0;%norm(Y0,'fro')/sqrt(n*m)*sqrt(0);
% Theta_corr = 0.2;
Test_num = 10;
Nl = [30,50,70,100];
pl = [0.5,1,1.5,2,2.5];
Error_l3_O = zeros(length(pl), length(Nl));
Error_l3_sphere = zeros(length(pl), length(Nl));

Theta = 0.5;
SNR = +inf;

%% Gen dict
for ni = 1:1:length(Nl)
    for mi = 1:1:length(pl)
        N = Nl(ni);
        K = N;
        M = round(10*N^pl(mi));
        Dini = randn(N);
        [Q,~] = qr(Dini);
        D0 = Q;
        Conv_l3_O_sum =zeros(1,80);
        Conv_l3_sph_sum =zeros(1,80);
        success_O = zeros(Test_num, 1);
        success_S = zeros(Test_num, 1);
        for it = 1:1:Test_num
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
            Y = Y0;
            %% ODL initial
            [dic_size, sample_size] = size(Y);
            [A0,~] = qr(randn(dic_size));
            A0 = A0(:,1:dic_size);
            %% L3 DL O
            A_l3_2stage = A0;
            last_A_l3_2stage = zeros(size(A0));
            e_l3_2stage = 1;
            iter_l3_2stage = 1;
            % f_val_l3_2stage = [];
            conv_l3_stage1 = [];
            while (e_l3_2stage>1e-4&&iter_l3_2stage < 100)
                AY_l3_2stage =  A_l3_2stage'*Y;
                dA_l3_2stage = (3 * abs(AY_l3_2stage) .*    AY_l3_2stage  * Y')';%AY.^(lp-1)*Y'/m;
                [U_l3_2stage,~,V_l3_2stage] = svd(-dA_l3_2stage,'econ');%-12*m*p^2*A%+2*p*noi_var^2+noi_var^4
                A_l3_2stage = U_l3_2stage*V_l3_2stage';
                [A_De]=De_permutation(A_l3_2stage,D0);
                Conv = norm( A_De  - D0,'fro')/norm(D0,'fro');
                %     f_val_l3_2stage = [f_val_l3_2stage sum(sum((abs(A_l3_2stage' * Y).^4)))/sample_size];
                e_l3_2stage = norm(A_De-last_A_l3_2stage,'fro')/(dic_size);
                conv_l3_stage1 = [conv_l3_stage1 Conv];
                last_A_l3_2stage = A_De;
                iter_l3_2stage =  iter_l3_2stage + 1;
                %    tl3 = toc
            end
            %stage 2
            r = last_A_l3_2stage;
            [n,p] = size(Y);
            conv_l3_stage2 = [];
            Obj = @(h) 1/(p)  *sum(sum(abs( h'*Y)));
            q_0 = r; %initialization
            q = q_0;
            mu = 1; rho = 0.9;  mu_threshold = 1e-10;
            iter = 1;
            while (iter < 51&&mu> mu_threshold)
                grad = 1/(p) *  Y*sign(Y'*q);
                
                q_new = r + 0.5*((q - mu*grad)-r*(q - mu*grad)'*r);
                obj_old = Obj(q); obj_new = Obj(q_new);
                grad_norm = norm(grad)^2;
                
                if obj_old <= obj_new && mu>mu_threshold
                    mu = mu*rho;
                elseif obj_old  >  obj_new && mu>mu_threshold
                    q = q_new; %update
                    iter = iter +1;
                    [U_l3_2stage,~,V_l3_2stage] = svd(q,'econ');%-12*m*p^2*A%+2*p*noi_var^2+noi_var^4
                    A_l3_2stage_2 = U_l3_2stage*V_l3_2stage';
                    [A_De_2]=De_permutation(A_l3_2stage_2,D0);
                    Conv_2 = norm( A_De_2  - D0,'fro')/norm(D0,'fro');
                    conv_l3_stage2 = [conv_l3_stage2   Conv_2];
                    %         if  nargin > 2
                    %             precond_q = real(ifft(  opts.precond .* fft(q)));
                    %             a_bar = normc(precond_q); a_bar = normc(real(ifft(1./fft(a_bar))));
                    %             dist2a = [dist2a dist_a(opts.a_0,a_bar)];
                    %         end
                    
                end
                
                
            end
            
            if conv_l3_stage2(end) <=1e-3   % found all standard basis
                success_O(it) = 1;
            end
            %         fprintf(' ODL proposed done!\n');
            
            %}
            %% sphere DL initial
            %}
            A_l3_2stage_sph = A0(:,1);
            last_A_l3_2stage_sph = zeros(size(  A_l3_2stage_sph));
            e_l3_2stage_sph = 1;
            iter_l3_2stage_sph = 1;
            % f_val_l3_2stage = [];
            conv_l3_stage1_sph = [];
            Conv_spht = 100;
            while ( e_l3_2stage_sph>1e-4&& iter_l3_2stage_sph < 100)
                AY_l3_2stage_sph =  A_l3_2stage_sph'*Y;
                dA_l3_2stage_sph = (3 * abs(AY_l3_2stage_sph) .* AY_l3_2stage_sph* Y')';%AY.^(lp-1)*Y'/m;
                %[U_l3_2stage_sph,~,V_l3_2stage] = svd(-dA_l3_2stage,'econ');%-12*m*p^2*A%+2*p*noi_var^2+noi_var^4
                A_l3_2stage_sph =dA_l3_2stage_sph/norm(dA_l3_2stage_sph);
                res_l3_sph = A_l3_2stage_sph'*D0;
                [~,ind]=max(abs(res_l3_sph));
                sign_t = sign(res_l3_sph(ind));
                Conv_sph = norm(A_l3_2stage_sph-sign_t*D0(:,ind),2);
                %     f_val_l3_2stage = [f_val_l3_2stage sum(sum((abs(A_l3_2stage' * Y).^4)))/sample_size];
                e_l3_2stage_sph = norm(A_l3_2stage_sph-last_A_l3_2stage_sph,2);
                % if Conv_sph> Conv_spht
                %     A_l3_2stage_sph = last_A_l3_2stage_sph;
                %end
                %res_l3_sph = A_l3_2stage_sph'*D0;
                %[~,ind]=max(abs(res_l3_sph));
                %sign_t = sign(res_l3_sph(ind));
                %Conv_sph = norm(A_l3_2stage_sph-sign_t*D0(:,ind),2);
                
                last_A_l3_2stage_sph = A_l3_2stage_sph;
                iter_l3_2stage_sph =  iter_l3_2stage_sph + 1;
                conv_l3_stage1_sph = [conv_l3_stage1_sph     Conv_sph ];
                Conv_spht = Conv_sph;
                %    tl3 = toc
            end
            
            r = last_A_l3_2stage_sph ;
            [n,p] = size(Y);
            conv_l3_stage2_sph= [];
            Obj = @(h) 1/(p)  *sum(sum(abs( h'*Y)));
            q_0 = r; %initialization
            q = q_0;
            mu = 1; rho = 0.85;  mu_threshold = 1e-10;
            iter = 1;
            while (iter < 51&&mu> mu_threshold)
                grad = 1/(p) *  Y*sign(Y'*q);
                mu = mu*rho;
                q_new = q - mu*grad;
                %projection
                q_new = q_new - (r'*q_new -1) /norm(q_new,'fro')^2 *r;
                obj_old = Obj(q); obj_new = Obj(q_new);
                grad_norm = norm(grad)^2;
                
                if obj_old <= obj_new && mu>mu_threshold
                    mu = mu*rho;
                elseif obj_old  >  obj_new && mu>mu_threshold
                    q = q_new; %update
                    iter = iter +1;
                    A_l3_2stage_2 =q/norm(q);
                    
                    %         if  nargin > 2
                    %             precond_q = real(ifft(  opts.precond .* fft(q)));
                    %             a_bar = normc(precond_q); a_bar = normc(real(ifft(1./fft(a_bar))));
                    %             dist2a = [dist2a dist_a(opts.a_0,a_bar)];
                    %         end
                    
                end
            end
            
            res_sph = A_l3_2stage_2'*D0;
            [~,ind]=max (abs(res_sph));
            sign_t = sign(res_sph(ind));
            Conv_2_sph = norm(A_l3_2stage_2 - sign_t*D0(:,ind),'fro')/norm( D0(:,ind),'fro');
            conv_l3_stage2_sph = [conv_l3_stage2_sph   Conv_2_sph ];
            
            if conv_l3_stage2_sph(end) <=1e-3   % found all standard basis
                success_S(it) = 1;
                
            end
            
            %                fprintf(' SDL proposed done!\n');
        end
        disp(['O: Success rate for N = ', num2str(N), 'M = ', num2str(pl(mi)), ': ', num2str(mean(success_O))]);
        disp(['S: Success rate for N = ', num2str(N), 'M = ', num2str(pl(mi)), ': ', num2str(mean(success_S))]);
        Error_l3_O(mi,ni) = mean(success_O);
        Error_l3_sphere(mi,ni) = mean(success_S);
    end
end

save('SAM_O_the05.mat','Error_l3_O');
save('SAM_SPH_the05.mat','Error_l3_sphere');

