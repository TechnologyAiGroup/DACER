function [flag,prediction_array,wModellist,ind_time,tot_time, try_num] = get_model_down(x_XPw,y_XPw,x,y,W,trainSetChallenges,predictionarray1,wModellist1,flag1,flag_diag,flag_earlystop,flag_fitloss,loss_func,loss_alpha,dis_mu,dis_sig)
%We simulate (x,y)-MXPUF
%MXPUF parameter
chalSize = 64;    % Bit length of challenge
mu = 0;           % Mean of variation in delay parameters
sigma = 1;        % Standard deviation of variation in delay parameters
        % x - number of APUFs in x-XOR PUF
        % y - number of APUFs in y-XOR PUF
feedback_a =33;   % feedback position to connect the output of x-XOR PUF and 
                  % the y-XOR PUF,  0<=feedback_a<=chalSize-1 


%We attack MXPUF U times and MXPUF is a noisy PUF with noise of sigmaNoise
%and check whether all APUFs x-XOR PUF can be modeled or not 
% x-XORPUF is the upper part, y-XOR PUF is the lower part

U=20;

%Time Evaluation and sigmaNoise
Evaluations =11;
sigmaNoise = 0.1;%noise_rate
nTrS = size(trainSetChallenges,1);
Size = chalSize+1;

%This array let us know how the occurence of found models matching some
%APUFs at x-XOR PUF: Xfound(1)->APUF(1), ..., Xfound(2) ->APUF(2).
wModellist=wModellist1;
prediction_array = predictionarray1;
flag = flag1;



    TrS= trainSetChallenges;
    challengeLR= TrS;
    challengePhiLR = Transform(challengeLR, nTrS, chalSize);
    challengePhiLR = fliplr(challengePhiLR);
    
    [Yp, ~] = classify(challengePhiLR,W,x);
    
    TrSp = zeros(nTrS,chalSize+1);
    for i=1:nTrS
        for j=1:(feedback_a-1)
            TrSp(i,j)= TrS(i,j);
        end
        TrSp(i,feedback_a)= Yp(i);
%         TrSp(i,feedback_a)= round(rand(1,1)*1);
        for j=(feedback_a+1):(chalSize+1)
            TrSp(i,j) = TrS(i,j-1);
        end                
    end
    
  Phi_TrSp = Transform(TrSp, nTrS, chalSize+1);
fl0=zeros(y,1);
fl=0;
time = zeros(U,1);

for k=1:U
    %fprintf('feedback_a %d-th and run %d-th \n',feedback_a, k);
    %Generate Traing Set, i.e., set of challenges
    [AResponse_TrS]= ComputeNoisyResponsesMXPUF(x_XPw,y_XPw,x,y,feedback_a,...
                                     TrS,size(TrS,1),chalSize,mu,sigma,sigmaNoise,Evaluations);
    AResponse_TrS = transpose(AResponse_TrS);
    InformReliability = ComputeTheNoiseInformation(AResponse_TrS,nTrS,Evaluations);
    
    
    %Since we focus on linear approximation attack at y-XOR PUF, we need to
    %modify the challenge
    %[~,wModel] = CMAES(Phi_TrSp,InformReliability,Size+1);
    opts=[];
    tic
    wModel = cmaes('modelAcc',rand(Size+1,1),0.5,opts, ...
            Phi_TrSp,InformReliability, ...
            flag_diag,flag_earlystop, ...
            flag_fitloss, loss_alpha, loss_func, dis_mu, dis_sig);
    time(k,1)=toc;
    % Compare wModel with APUFs at y-XOR PUF. 
        
    result=zeros(1,y); 
  
    resultpre=zeros(1,y);
    for i=1:y
         wp = zeros(1,Size+1);
         for j=1:(Size+1)
             wp(j)=y_XPw(i,j);
         end            
         result(i)=CompareTwoModels(wp,wModel,Size+1);
         if(result(i)>0.9||result(i)<0.1)
             fl0(i,1)=1;
         end
         if(result(i)<1-result(i))
             resultpre(i)=1-result(i);
         else
             resultpre(i)=result(i);
         end
    end 
%     [~,index]=max(resultpre);
    for index=1:y
        if(result(index)~=resultpre(index))
        if(prediction_array(index)<resultpre(index))
            prediction_array(index)=resultpre(index);
            flag(index)=1;
            wModel=transpose(wModel);
            for j=1:(Size+1)
              wModellist(index,j)=-wModel(j);
            end  
        end
      else
         if(prediction_array(index)<resultpre(index))
            prediction_array(index)=resultpre(index);
            flag(index)=0;
            wModel=transpose(wModel);
            for j=1:(Size+1)
              wModellist(index,j)=wModel(j);
            end  
         end
        
       end
    end
    for p=1:y
        if(fl0(p,1)==0)
            break;
        end
        if(p==y)
            fl=1;
        end
    end
    if fl==1
        break;
    end
end  
       ind_time = sum(time) / k;
       tot_time = sum(time) / 60;
       try_num = k;
       disp(['Ind time = ', num2str(ind_time), ...
           ' s, Tot time = ', num2str(tot_time), ' min, Try num = ', num2str(try_num)])
end
