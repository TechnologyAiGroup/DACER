function [NoiseInform] = ComputeTheNoiseInformation(AResponses,nRows,Evaluations)

%We are given an array of XORPUF after Evaluations measurement, i.e.,
%AResponses[i,j] is the output of XORPUF of challenge[i] and Evaluation[j].
%We compute the information of the reliability NoiseInform[i]=|Evaluations/2 - \sum AResponse[i,e]| 

NoiseInform = zeros(1,nRows);

for i=1:nRows
  
   for e=1:Evaluations
       NoiseInform(i) = NoiseInform(i) + AResponses(i,e);
   end
   
   NoiseInform(i) = abs((Evaluations/2)-NoiseInform(i));
   %NoiseInform(i) = abs(NoiseInform(i)-5)*2; 
   %NoiseInform(i) = NoiseInform(i)/10;
     
end

end

