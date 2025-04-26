function [AResponse,select] = ComputeResponsecMPUF(Dw,Sw,nXOR,APhi,nRows,Size)
% The function computes the array of responses for a given array of feature
% vectors APhi and weight vector w
% The response = 0 if  w*APhi >0, otherwise = 1;
%   Detailed explanation goes here
  
  AResponse = zeros(1,nRows);
  select = zeros(1,nRows);
  for i=1:nRows
      %Outputs of each APUF instances 
      Sum = zeros(1,nXOR);
      
      %Compute the outputs of nXOR APUFs
      for j=1:Size
          for k = 1:nXOR
            Sum(k) = Sum(k) + Sw(k,j)*APhi(i,j);  
          end
      end     
      for x=1:nXOR
          if Sum(x)>=0
             Sum(x)=1;
          else
             Sum(x)=0;
          end
          select(i)=select(i)*2+Sum(x);
          
      end
      %Compute the output of XORPUF 
      temp=(select(i)+2)/2;
      temp=floor(temp);
      neg=mod((select(i)+2),2);
      for j=1:Size
            AResponse(i) = AResponse(i) + Dw(temp,j)*APhi(i,j);  
      end  
      if(AResponse(i)>=0)
          AResponse(i)=1;
      else 
          AResponse(i)=0;
      end
      if(neg==1)
          AResponse(i)=~AResponse(i);
      end
  end

end

