%We attack MXPUF based on Becker's attack for one PUF at y-XOR PUF and
%y-XOR PUF. We do the linear approximation attack. 
function [flag,prediction_array,wModellist] = MXPUF_getfirstmodel(x_XPw,y_XPw,x,y,trainSetChallenges,flag1,prediction_array1,wModellist1,flag_diag,flag_earlystop,flag_fitloss,loss_func,loss_alpha,dis_mu,dis_sig)
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

% U=35;
U=20;
%Time Evaluation and sigmaNoise
Evaluations =11;
sigmaNoise = 0.1;%noise_rate
nTrS = size(trainSetChallenges,1);
Size = chalSize+1;
%earlystop=0;
%This array let us know how the occurence of found models matching some
%APUFs at x-XOR PUF: Xfound(1)->APUF(1), ..., Xfound(2) ->APUF(2).

prediction_array = prediction_array1;
flag = flag1;
wModellist=wModellist1;
TrS= trainSetChallenges;
TrSp = zeros(nTrS,chalSize+1);
        for i=1:nTrS
            for j=1:(feedback_a-1)
                TrSp(i,j)= TrS(i,j);
            end
            randomshu= round(rand(1,1)*1);
            TrSp(i,feedback_a)= randomshu;
            for j=(feedback_a+1):(chalSize+1)
                TrSp(i,j) = TrS(i,j-1);
            end                
        end

  Phi_TrSp = Transform(TrSp, nTrS, chalSize+1);
  if(feedback_a==33)  
  %if((feedback_a==1)||(feedback_a==chalSize/2)||(feedback_a==chalSize))   
    for k=1:U
    %fprintf('feedback_a %d-th and run %d-th \n',feedback_a, k);
    %Generate Traing Set, i.e., set of challenges
    [AResponse_TrS]= ComputeNoisyResponsesMXPUF(x_XPw,y_XPw,x,y,feedback_a,...
                                     TrS,nTrS,chalSize,mu,sigma,sigmaNoise,Evaluations);
    AResponse_TrS = transpose(AResponse_TrS);
    InformReliability = ComputeTheNoiseInformation(AResponse_TrS,nTrS,Evaluations);
    
    
    %Since we focus on linear approximation attack at y-XOR PUF, we need to
    %modify the challenge
    
    %[~,wModel] = CMAES(Phi_TrSp,InformReliability,Size+1);
    opts=[];
    wModel = cmaes('modelAcc',rand(Size+1,1),0.5,opts, ...
            Phi_TrSp,InformReliability, ...
            flag_diag,flag_earlystop, ...
            flag_fitloss, loss_alpha, loss_func, dis_mu, dis_sig);
    % Compare wModel with APUFs at y-XOR PUF. 
    result=zeros(1,y); 
    resultpre=zeros(1,y);
    for i=1:y
         wp = zeros(1,Size+1);
         for j=1:(Size+1)
             wp(j)=y_XPw(i,j);
         end            
         result(i)=CompareTwoModels(wp,wModel,Size+1);
         if(result(i)<1-result(i))
             resultpre(i)=1-result(i);
         else
             resultpre(i)=result(i);
         end
    end 
   % [~,index]=max(resultpre);
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
    countacc=0;
    for i=1:y
       if prediction_array(i)>=0.95
           countacc=countacc+1;
       end
    end
    if countacc ==y
        break;
    end
    
    end 
  end   
    
end  

    




