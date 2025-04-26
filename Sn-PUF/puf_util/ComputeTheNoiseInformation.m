function [NoiseInform] = ComputeTheNoiseInformation(AResponses,nRows,Evaluations)

NoiseInform = zeros(1,nRows);

for i=1:nRows
  
   for e=1:Evaluations
       NoiseInform(i) = NoiseInform(i) + AResponses(i,e);
   end
   
   NoiseInform(i) = abs((Evaluations/2)-NoiseInform(i));
     
end

end

