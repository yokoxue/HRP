function[A_Spams,tot_Spams] = Dic_Learn_Spams (Y,A0)

    %% DL
  
    param.D =  A0;
    param.lambda=0.15;
    param.numThreads=-1; % number of threads
    param.batchsize=400;
    param.verbose=false;
    param.modeD=0;

    param.iter=1000;  % let us see what happens after 1000 iterations.
    tic
    A_Spams = mexTrainDL(Y,param);
    tot_Spams=toc;
    
   
