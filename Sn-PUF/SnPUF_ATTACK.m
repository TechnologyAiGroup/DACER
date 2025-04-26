function [acc,tot_time] = SnPUF_ATTACK(chalSize,nXOR,nTrS,nTeS,sigmaNoise, ...
                    flag_diag,flag_earlystop, ...
                    flag_fitloss,loss_func,loss_alpha,dis_mu,dis_sig)

mu = 0;           % Mean of variation in delay parametersd
sigma = 1;     % Standard deviation of variation in delay parameters  

Size = chalSize+1;

Evaluations =11;

XORw= XORPUFgeneration(nXOR*2,chalSize,mu,sigma);
randa= randperm(nXOR/2,nXOR/2)+(nXOR/2);


TeS= randi([0 1], nTeS, chalSize);
TeS1(:,1:chalSize/2)= TeS(:,(chalSize/2)+1:chalSize);
TeS1(:,(chalSize/2)+1:chalSize)= TeS(:,1:chalSize/2);
Phi_TeS = Transform(TeS, nTeS, chalSize);
Phi_TeS1 = Transform(TeS1, nTeS, chalSize);
for i=1:nXOR
    AResponse_TeS(i,:) = ComputeResponseXOR(XORw((i-1)*2+1,:),1,Phi_TeS,nTeS,Size);
    AResponse_TeS1(i,:) = ComputeResponseXOR(XORw((i-1)*2+2,:),1,Phi_TeS1,nTeS,Size);
    AResponse_TeS(i,:) =xor(AResponse_TeS(i,:),AResponse_TeS1(i,:));
end
for i=1:(nXOR/2)
    Response_TeS(i,:)=and(AResponse_TeS(i,:),AResponse_TeS(randa(i),:));
end
Response_TeS1(1,:)=Response_TeS(1,:);
for i=2:(nXOR/2)
    Response_TeS1(1,:)=xor(Response_TeS1(1,:),Response_TeS(i,:));
end

TrS= randi([0 1], nTrS, chalSize);
TrS1(:,1:chalSize/2)= TrS(:,(chalSize/2)+1:chalSize);
TrS1(:,(chalSize/2)+1:chalSize)= TrS(:,1:chalSize/2);
Phi_TrS = Transform(TrS, nTrS, chalSize);
Phi_TrS1 = Transform(TrS1, nTrS, chalSize);
for p=1:Evaluations
for i=1:nXOR
    AResponse_TrS(i,:) = ComputeNoisyResponsesXOR(XORw((i-1)*2+1,:),1,Phi_TrS,nTrS,Size,chalSize,sigma,sigmaNoise,1);
    AResponse_TrS1(i,:) = ComputeNoisyResponsesXOR(XORw((i-1)*2+2,:),1,Phi_TrS,nTrS,Size,chalSize,sigma,sigmaNoise,1);
    AResponse_TrS(i,:) =xor(AResponse_TrS(i,:),AResponse_TrS1(i,:));
end
for i=1:(nXOR/2)
    Response_TrS(i,:)=and(AResponse_TrS(i,:),AResponse_TrS(randa(i),:));
end
Response_TrS1(1,:)=Response_TrS(1,:);
for i=2:(nXOR/2)
    Response_TrS1(1,:)=xor(Response_TrS1(1,:),Response_TrS(i,:));
end
AResponse_TrS_final(p,:)=Response_TrS1;
end

model=zeros(nXOR,Size);
InformReliability = ComputeTheNoiseInformation(AResponse_TrS_final',nTrS,Evaluations);
%cma-es attack
opts=[];
tic
wModel = cmaes('modelAcc',rand((Size+1)*nXOR,1),0.5,opts,Phi_TrS,InformReliability, flag_diag,flag_earlystop, ...
            flag_fitloss, loss_alpha, loss_func, dis_mu, dis_sig,AResponse_TrS,nXOR,Size,AResponse_TrS_final,randa,Phi_TrS1);
tot_time=toc;

for i=1:nXOR
    model(i,:)=wModel((i-1)*(Size+1)+1:(i-1)*(Size+1) + Size,:);
end
acc=modelAcc(model,Phi_TeS,Response_TeS1,Size,nTeS,nXOR);
acc = max(acc, 1 - acc);
