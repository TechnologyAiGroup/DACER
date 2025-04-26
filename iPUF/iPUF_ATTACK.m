function [ac,ind_time,tot_time, try_num] = iPUF_ATTACK(chalSize,nhighXOR,nbelowXOR,nTrS,nTest,sigmaNoise, ...
                    flag_diag,flag_earlystop, ...
                    flag_fitloss,loss_func,loss_alpha,dis_mu,dis_sig)
mu = 0;           % Mean of variation in delay parameters
sigma = 1;        % Standard deviation of variation in delay parameters
x =nhighXOR;        % x - number of APUFs in x-XOR PUF
y =nbelowXOR;        % y - number of APUFs in y-XOR PUF
feedback_a =33;   % feedback position to connect the output of x-XOR PUF and 
                  % the y-XOR PUF,  0<=feedback_a<=chalSize-1 
Size = chalSize+1;

%generate (x,y)-MXPUF
prediction_array = zeros(1,y);
flag = zeros(1,y);
wModellist=zeros(y,Size+1);
[x_XPw,y_XPw]=MXPUFgeneration(x,y,chalSize,mu,sigma);

%generate Test Set and Training Set. 
xunhuancount=1;
zongxunhuan=0;
trainSetChallenges= randi([0 1], nTrS, chalSize);
trainSetResponses = ComputeResponseMXPUF( ...
                       x_XPw,y_XPw,x,y,feedback_a, ...
                       trainSetChallenges,nTrS,chalSize ...
                       );  

% x_XPw=csvread("./dataset/1_3_64_0.10_Xw.csv");
% y_XPw=csvread("./dataset/1_3_64_0.10_Yw.csv");
% trainSetChallenges=csvread("./dataset/1_3_64_0.10_Challenge.csv");
% trainSetResponses=csvread("./dataset/1_3_64_0.10_Response.csv");

[flag,prediction_array,wModellist] = MXPUF_getfirstmodel(x_XPw,y_XPw,x,y,trainSetChallenges,flag,prediction_array,wModellist,flag_diag,flag_earlystop,flag_fitloss,loss_func,loss_alpha,dis_mu,dis_sig);

[testchallenge,testresponse]=gettest_set(x_XPw,chalSize,x);
preac=0;
preWup=zeros(x,chalSize+1);
while xunhuancount<4&&zongxunhuan<24

zongxunhuan=zongxunhuan+1;

if(preac>0.8)
    xunhuancount=xunhuancount+1;
end

TrSp = zeros(nTrS,chalSize+1);
    for i=1:nTrS
        for j=1:(feedback_a-1)
            TrSp(i,j)= trainSetChallenges(i,j);
        end
        Interposbit= 0;
        TrSp(i,feedback_a)= Interposbit;
        for j=(feedback_a+1):(chalSize+1)
            TrSp(i,j) = trainSetChallenges(i,j-1);
        end                
    end
    Phi_TrSp = Transform(TrSp, nTrS, chalSize+1);
    
   Areponse0=ComputeResponseXOR(wModellist,y,Phi_TrSp,nTrS,size(Phi_TrSp,2));
   
TrSp = zeros(nTrS,chalSize+1);
    for i=1:nTrS
        for j=1:(feedback_a-1)
            TrSp(i,j)= trainSetChallenges(i,j);
        end
        Interposbit= 1;
        TrSp(i,feedback_a)= Interposbit;
        for j=(feedback_a+1):(chalSize+1)
            TrSp(i,j) = trainSetChallenges(i,j-1);
        end                
    end
    Phi_TrSp = Transform(TrSp, nTrS, chalSize+1);
    Areponse1=ComputeResponseXOR(wModellist,y,Phi_TrSp,nTrS,size(Phi_TrSp,2));
    LRtrainsetchallenges=[];
    LRtrainsetresponse=[];
    count=1;
    for i=1:nTrS
        if(Areponse0(i)~= Areponse1(i))
            if(Areponse0(i)==trainSetResponses(i))
                LRtrainsetchallenges(count,:)=trainSetChallenges(i,:);
                LRtrainsetresponse=[LRtrainsetresponse,0];
            else
                LRtrainsetchallenges(count,:)=trainSetChallenges(i,:);
                LRtrainsetresponse=[LRtrainsetresponse,1];
            end
            count=count+1;
        end
    end
LRtrainsetresponse=transpose(LRtrainsetresponse);
[Wup,ac1]= get_model_up(LRtrainsetchallenges,LRtrainsetresponse,x,testchallenge,testresponse,preac,preWup);

preac=ac1;
preWup=Wup;
[flag,prediction_array,wModellist,ind_time,tot_time, try_num] = get_model_down(x_XPw,y_XPw,x,y,Wup,trainSetChallenges,prediction_array,wModellist,flag,flag_diag,flag_earlystop,flag_fitloss,loss_func,loss_alpha,dis_mu,dis_sig);

end

TestSetChallenges= randi([0 1], nTest, chalSize);
TestSetResponses = ComputeResponseMXPUF( ...
                       x_XPw,y_XPw,x,y,feedback_a, ...
                       TestSetChallenges,nTest,chalSize ...
                       ); 

TestY = ComputeResponseMXPUF2( ...
                       Wup,wModellist,x,y,feedback_a, ...
                       TestSetChallenges,nTest,chalSize ...
                       ); 
                   
[ac, precision, recall, fscore] = accuracy(TestSetResponses,TestY);



