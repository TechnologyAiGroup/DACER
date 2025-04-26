function [acc,ind_time,tot_time, try_num] = MPUF_ATTACK(chalSize,nXOR,nTrS,nTeS,sigmaNoise, ...
                    flag_diag,flag_earlystop, ...
                    flag_fitloss,loss_func,loss_alpha,dis_mu,dis_sig)
    acc = 0.0;
    ind_time = 0.0;
    tot_time = 0.0;
    try_num = 0.0;
    
    mu = 0;           % Mean of variation in delay parametersd
    sigma = 1;     % Standard deviation of variation in delay parameters
    Sw= XORPUFgeneration(nXOR,chalSize,mu,sigma);
    Dw= XORPUFgeneration(2^nXOR,chalSize,mu,sigma);
    
    zong=nTrS+nTeS;
    Size = chalSize+1;
    
    Evaluations =11;
    
    U=50;
    matchingRateMultipleTimes = zeros(U,nXOR);
    
    zongS= randi([0 1], zong, chalSize);
    Phi_zongS = Transform(zongS, zong, chalSize);
    AResponse_zongS = ComputeNoisyResponsesMPUF(Dw,Sw,nXOR,Phi_zongS,zong,Size,chalSize,sigma,sigmaNoise,Evaluations);
    
    % Phi_zongS=csvread("./dataset/3_64_0.10_Phi.csv");
    % AResponse_zongS=csvread("./dataset/3_64_0.10_Response.csv");
    % Sw=csvread("./dataset/3_64_0.10_Sw.csv");

    Phi_TrS = Phi_zongS(1:nTrS,:);
    AResponse_TrS=AResponse_zongS(1:nTrS,:);
    Phi_TeS = Phi_zongS(nTrS+1:zong,:);
    AResponse_TeS=AResponse_zongS(nTrS+1:zong,:);
    
    for i=1:nTeS
        ab=0;
        for j=1:Evaluations
            ab=ab+AResponse_TeS(i,j);
        end
        if ab>5.5
            AResponse_TeS2(1,i)=1;
        else
            AResponse_TeS2(1,i)=0;
        end
    end
    
    modelS=zeros(nXOR,Size);
    modelD=zeros(2^nXOR,Size);
    time = zeros(U,1);
    num=zeros(2^nXOR,1);
    CD=zeros(nXOR,nTrS,Size);
    RD=zeros(nXOR,nTrS,1);
    AResponseSAPUFs = zeros(nXOR,nTeS);
    for i=1:nXOR 
        AResponseSAPUFs(i,:)=ComputeResponseXOR(Sw(i,:),1,Phi_TeS,nTeS,Size);
    end
    
    InformReliability = ComputeTheNoiseInformation(AResponse_TrS,nTrS,Evaluations);
    
    flag=0;
    modelacc=zeros(1,2^nXOR);   
    % CMA-ES ATTACK find the Selection APUF
    for k=1:U
        opts=[];
        tic
        wModel = cmaes('modelAcc',rand(Size,1),0.5,opts, ...
                Phi_TrS,InformReliability, ...
                flag_diag,flag_earlystop, ...
                flag_fitloss, loss_alpha, loss_func, dis_mu, dis_sig);
        time(k,1)=toc;
        for i=1:nXOR
           matchingRateMultipleTimes(k,i)=modelAccHa(wModel,Phi_TeS,AResponseSAPUFs(i,:),Size,nTeS);
           if((matchingRateMultipleTimes(k,i)>0.9)||(matchingRateMultipleTimes(k,i)<0.1))
               modelS(i,:)=wModel;
           end
        end

        for i=1:nXOR
            if modelS(i,1)==0
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
    % LR ATTACK find the Data APUF
    if flag==1
        select=ComputeselectMPUF(modelS,nXOR,Phi_TrS,nTrS,Size);
        for i=1:nTrS
            temp=select(i)+1;
            num(temp)=num(temp)+1;
            temp1=num(temp);
            CD(temp,temp1,:)=Phi_TrS(i,:);
            RD(temp,temp1,:)=AResponse_TrS(i,1);  
        end
        for i=1:2^nXOR
            temp=num(i);
            tCD=CD(i,1:temp,:);
            tCD=reshape(tCD,temp,Size);
            tRD=RD(i,1:temp);
            tRD=tRD';
            split=floor(temp*0.9);
            [allac, precision, recall, fscore,W]=LR_XAPUF(tCD(1:split,:),tRD(1:split,:),tCD(split+1:end,:),tRD(split+1:end,:),chalSize,1);
            modelD(i,:)=W;
    
        end
        [AResponse_TeS1,select] = ComputeResponseMPUF(modelD,modelS,nXOR,Phi_TeS,nTeS,Size);
        acc=0;
        for i=1:nTeS
            if AResponse_TeS1(i)==AResponse_TeS2(i)
                acc=acc+1;
            end
        end
        acc=acc/nTeS;
        acc=max(acc,1-acc);
        ind_time = sum(time) / k;
        tot_time = sum(time) / 60;
        try_num = k;
        disp(['Accuracy = ', num2str(acc), ', Ind time = ', num2str(ind_time), ...
               ' s, Tot time = ', num2str(tot_time), ' min, Try num = ', num2str(try_num)])
    else
          disp('Failed') 
    end
end


