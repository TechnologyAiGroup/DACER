function [acc,ind_time,tot_time, try_num] = XORPUF_ATTACK(chalSize,nXOR,nTrS,nTeS,sigmaNoise, ...
                    flag_diag,flag_earlystop, ...
                    flag_fitloss,loss_func,loss_alpha,dis_mu,dis_sig)
    acc = 0.0;
    ind_time = 0.0;
    tot_time = 0.0;
    try_num = 0.0;

    mu = 0;           % Mean of variation in delay parametersd
    sigma = 1;     % Standard deviation of variation in delay parameters  
    
    Size = chalSize + 1;
    
    Evaluations = 11;
    
    U = 50;
    
    matchingRateMultipleTimes = zeros(U,nXOR);
    XORw = XORPUFgeneration(nXOR,chalSize,mu,sigma);
    
    TeS= randi([0 1], nTeS, chalSize);
    Phi_TeS = Transform(TeS, nTeS, chalSize);
    AResponse_TeS1 = ComputeResponseXOR(XORw,nXOR,Phi_TeS,nTeS,Size);
    
    TrS= randi([0 1], nTrS, chalSize);
    Phi_TrS = Transform(TrS, nTrS, chalSize);
    AResponse_TrS = ComputeNoisyResponsesXOR(XORw,nXOR,Phi_TrS,nTrS,Size,chalSize,sigma,sigmaNoise,Evaluations);
    InformReliability = ComputeTheNoiseInformation(AResponse_TrS,nTrS,Evaluations);
    
    AResponseALLAPUFs = zeros(nXOR,nTeS);
    for i=1:nXOR 
        AResponseALLAPUFs(i,:) = ComputeResponseXOR(XORw(i,:),1,Phi_TeS,nTeS,Size);
    end
    %
    %cma-es attack
    flag = 0;
    time = zeros(U,1);
    model = zeros(nXOR,Size);

    % Phi_TrS=csvread("./dataset/4_64_0.10_Phi_TrS.csv");
    % InformReliability=csvread("./dataset/4_64_0.10_InformReliability.csv");
    % Phi_TeS=csvread("./dataset/4_64_0.10_Phi_TeS.csv");
    % AResponse_TeS1=csvread("./dataset/4_64_0.10_TeS_Response.csv");
    % AResponseALLAPUFs=csvread("./dataset/4_64_0.10_All_Response.csv");
    
    for k = 1:U
        opts=[];
        tic
        wModel = cmaes('modelAcc',rand(Size,1),0.5,opts, ...
            Phi_TrS,InformReliability, ...
            flag_diag,flag_earlystop, ...
            flag_fitloss, loss_alpha, loss_func, dis_mu, dis_sig);
        time(k,1)=toc;
        for i=1:nXOR
           matchingRateMultipleTimes(k,i)=modelAccHa(wModel,Phi_TeS,AResponseALLAPUFs(i,:),Size,nTeS);
           if((matchingRateMultipleTimes(k,i)>0.9)||(matchingRateMultipleTimes(k,i)<0.1))
               model(i,:)=wModel;
           end
        end
        for i=1:nXOR
            if model(i,1)==0
                break;
            end
            if i==nXOR
                flag=1;
            end
        end
        if flag ==1
           break;
        end
        
    end 
    if(flag==1)
       acc = modelAcc(model,Phi_TeS,AResponse_TeS1,Size,nTeS,nXOR);
       acc = max(acc, 1 - acc);
       ind_time = sum(time) / k;
       tot_time = sum(time) / 60;
       try_num = k;
       disp(['Accuracy = ', num2str(acc), ', Ind time = ', num2str(ind_time), ...
           ' s, Tot time = ', num2str(tot_time), ' min, Try num = ', num2str(try_num)])
    else
        disp('Failed')
    end
end