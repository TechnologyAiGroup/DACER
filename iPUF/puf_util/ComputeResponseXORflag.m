function [AResponse] = ComputeResponseXORflag(XORw,nXOR,APhi,nRows,Size,flag)
% The function computes the array of responses for a given array of feature
% vectors APhi and weight vector w
% The response = 0 if  w*APhi >0, otherwise = 1;
%   Detailed explanation goes here
  
  AResponse = ones(1,nRows);

  for i=1:nRows
      %Outputs of each APUF instances 
      Sum = zeros(1,nXOR);
      
      %Compute the outputs of nXOR APUFs
      for j=1:Size
          for k = 1:nXOR
            Sum(k) = Sum(k) + XORw(k,j)*APhi(i,j);  
          end
      end     
      
      %Compute the output of XORPUF 
      Prod = 1;
      for jk=1:nXOR
          if(flag(i)==1)
              Sum(i)=-Sum(i);
          end
      end
      for k=1:nXOR
          Prod = Prod*Sum(k);
      end

      if(Prod>0)
          AResponse(i)=0;
      else 
          AResponse(i)=1;
      end
  end

end