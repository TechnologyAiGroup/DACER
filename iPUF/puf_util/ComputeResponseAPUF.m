function [AResponse] = ComputeResponseAPUF(w,APhi,nRows,Size)
% The function computes the array of responses for a given array of feature
% vectors APhi and weight vector w
% The response = 0 if  w*APhi >0, otherwise = 1;
%   Detailed explanation goes here
  
  AResponse = ones(1,nRows);

  for i=1:nRows
      %Outputs of APUF instances 
      Sum = 0;
      
      for j=1:Size
          Sum = Sum + w(j)*APhi(i,j);
      end
      
      if(Sum>0)
          AResponse(i)=0;
      else 
          AResponse(i)=1;
      end
  end

end

