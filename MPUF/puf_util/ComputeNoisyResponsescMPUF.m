function [select,AResponse] = ComputeNoisyResponsescMPUF(Dw,Sw,nXOR,APhi,nRows,Size,ChalSize,sigma,sigmaNoise,Evaluations)


%We are given the APUF which has: XORw weight vectors, nXOR APUFs, APhi
%vectors computed from array of challenges and the array has nRows, Size
%columns. The APUF is created by a given sigma, mu and has ChalSize-bit
%challenge

%We need to create a noisy CRPs with sigmaNoise and for each Challenge, we
%evaluate it Evaluation times. It means that we have an array of responses
%AResponse having nRows and Evaluations column. Typically, AResponse[i,j]
%is the response of challange[i] of Evaluation[j]. 



% The function computes the array of responses for a given array of feature
% vectors APhi and weight vector w
% The response = 0 if  w*APhi >0, otherwise = 1;
%   Detailed explanation goes here
  
%Creating an array of noisy Responses AResponse with nRows rows and
%Evaluations columns. 

  AResponse = zeros(nRows,Evaluations);
  select = zeros(nRows,Evaluations);
%
  for e=1:Evaluations
      
      %Create a noise vector 
      %sigNoise = sigmaNoise*sigma;
      SXORwNoise= XORPUFgeneration(nXOR,ChalSize,0,sigmaNoise);
      %DXORwNoise=XORPUFgeneration(2^nXOR,ChalSize,0,sigNoise);
      %Create a XOR affected by the noise newXORw
      newSXORw = zeros(nXOR,Size);
      %newDXORw = zeros(2^nXOR,Size);
      for j=1:Size
          for k=1:nXOR
              newSXORw(k,j) = Sw(k,j) + SXORwNoise(k,j);
          end
      end
      %Evaluate the responses
      
      for i=1:nRows
           %Outputs of each APUF instances 
            Sum = zeros(1,nXOR);
      
           %Compute the outputs of nXOR APUFs
      
            for j=1:Size
               for k = 1:nXOR
                  Sum(k) = Sum(k) + newSXORw(k,j)*APhi(i,j);  
               end
            end     
            for x=1:nXOR
                if Sum(x)>=0
                    Sum(x)=1;
                else
                    Sum(x)=0;
                end
                select(i,e)=select(i,e)*2+Sum(x);
            end
            temp=(select(i,e)+2)/2;
            temp=floor(temp);
            neg=mod((select(i,e)+2),2);
           %Compute the output of XORPUF 
            for j=1:Size
                AResponse(i,e) = AResponse(i,e) + (Dw(temp,j)+normrnd(0,sigmaNoise,1,1))*APhi(i,j);  
            end  
            if(AResponse(i,e)>=0)
               AResponse(i,e)=1;
            else 
               AResponse(i,e)=0;
            end
            if(neg==1)
                AResponse(i,e)=~AResponse(i,e);
            end
       end
      
  end
  
   

end

