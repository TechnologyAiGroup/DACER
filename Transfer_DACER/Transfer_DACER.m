clear all;
clc;
addpath(genpath('./puf_util'));
chalSize = 64;    % Bit length of challenge
mu = 0;           % Mean of variation in delay parametersd
sigma = 1;     % Standard deviation of variation in delay parameters  

nXOR = 2;

nTrS =40000;
nTeS =500;

nTrS_t=4000;

Size = chalSize+1;


Evaluations =11;
sigmaNoise = 0.1;
U=50;

matchingRateMultipleTimes = zeros(U,nXOR);
matchingRateMultipleTimes_t = zeros(U,nXOR);
XORw= XORPUFgeneration(nXOR,chalSize,mu,sigma);
XORw_t= XORPUFgeneration(nXOR,chalSize,mu,sigma);

TrS= randi([0 1], nTrS, chalSize);
Phi_TrS = Transform(TrS, nTrS, chalSize);
AResponse_TrS = ComputeNoisyResponsesXOR(XORw,nXOR,Phi_TrS,nTrS,Size,chalSize,sigma,sigmaNoise,Evaluations);
AResponse_TrS_t = ComputeNoisyResponsesXOR(XORw_t,nXOR,Phi_TrS,nTrS,Size,chalSize,sigma,sigmaNoise,Evaluations);
%PCA
[coeff, Phi_TrS1, ~, ~, explained, mu] = pca(Phi_TrS);
Phi_TrS = [Phi_TrS Phi_TrS1];


InformReliability = ComputeTheNoiseInformation(AResponse_TrS,nTrS,Evaluations);
InformReliability_t = ComputeTheNoiseInformation(AResponse_TrS_t,nTrS,Evaluations);


%to verify the found APUF
AResponseALLAPUFs = zeros(nXOR,nTeS);
AResponseALLAPUFs_t = zeros(nXOR,nTeS);


TeS= randi([0 1], nTeS, chalSize);
Phi_TeS = Transform(TeS, nTeS, chalSize);
AResponse_TeS1 = ComputeResponseXOR(XORw,nXOR,Phi_TeS,nTeS,Size);
AResponse_TeS1_t = ComputeResponseXOR(XORw_t,nXOR,Phi_TeS,nTeS,Size);
X_test_centered = Phi_TeS - mu;  
%PCA
Phi_TeS1 = X_test_centered * coeff;  
Phi_TeS =  [Phi_TeS Phi_TeS1];

flag=0;
model=zeros(nXOR,Size*2);
C_all=zeros(nXOR,Size*2,Size*2);
xmean_all=zeros(nXOR,Size*2);
model_t=zeros(nXOR,Size*2);

for i=1:nXOR 
    AResponseALLAPUFs(i,:)=ComputeResponseXOR(XORw(i,:),1,Phi_TeS,nTeS,Size);
end
for i=1:nXOR 
    AResponseALLAPUFs_t(i,:)=ComputeResponseXOR(XORw_t(i,:),1,Phi_TeS,nTeS,Size);
end

% Phi_TrS=csvread("./dataset/2_64_0.10_Phi_TrS.csv");
% Phi_TeS=csvread("./dataset/2_64_0.10_Phi_TeS.csv");
% InformReliability=csvread("./dataset/2_64_0.10_Inform.csv");
% InformReliability_t = csvread("./dataset/2_64_0.10_Inform_t.csv");
% AResponse_TeS1_t=csvread("./dataset/2_64_0.10_Response_TeS.csv");
% AResponseALLAPUFs=csvread("./dataset/2_64_0.10_AllPUFs.csv");
% AResponseALLAPUFs_t=csvread("./dataset/2_64_0.10_AllPUFs_t.csv");

%cma-es attack for source puf
for k=1:U
    opts=[];
    [arfitness,arx,C,xmean] = CMAES(Phi_TrS,InformReliability,size(Phi_TrS,2));
    for i=1:nXOR
       matchingRateMultipleTimes(k,i)=modelAccHa(arx(:,1),Phi_TeS,AResponseALLAPUFs(i,:),size(Phi_TeS,2),nTeS);
       if((matchingRateMultipleTimes(k,i)>0.9)||(matchingRateMultipleTimes(k,i)<0.1))
           model(i,:)=arx(:,1);
           C_all(i,:,:)=C;
           xmean_all(i,:)=xmean;
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
    if flag == 1
       break;
    end
end 
%transfer
if flag == 1
flag = 0;
for k=1:U
    r=randi(nXOR);
    C_t=squeeze(C_all(r,:,:));
    mean_t=xmean_all(r,:)';
    [arfitness1,arx1,C1,mean1] = WS_CMAES(Phi_TrS(1:nTrS_t,:),InformReliability_t(:,1:nTrS_t),size(Phi_TrS,2),mean_t,0.1,C_t);
    for i=1:nXOR
       matchingRateMultipleTimes_t(k,i)=modelAccHa(arx1(:,1),Phi_TeS,AResponseALLAPUFs_t(i,:),size(Phi_TeS,2),nTeS);
       if((matchingRateMultipleTimes_t(k,i)>0.9)||(matchingRateMultipleTimes_t(k,i)<0.1))
           model_t(i,:)=arx1(:,1);
       end
    end
    for i=1:nXOR
        if model_t(i,1)==0
            break;
        end
        if i==nXOR
            flag=1;
        end
    end
    if flag == 1
       break;
    end
end
if flag == 1
    acc = modelAcc(model_t,Phi_TeS,AResponse_TeS1_t,Size*2,nTeS,nXOR);
    acc = max(acc, 1 - acc);
end
end