function [Rate] = ComputeTheNoise(AResponses,nRows,Evaluations,Threshold)

%We are given an array of XORPUF after Evaluations measurement, i.e.,
%AResponses[i,j] is the output of XORPUF of challenge[i] and Evaluation[j].
%We compute how many unreliable outputs wrt Threshold. 

Rate = 0;

for i=1:nRows
   %create a counters for counting 0s and 1s in Evaluations measurements
   count = zeros(1,2);
   %count(1) for 0 and count(2) for 1
   
   for e=1:Evaluations
      if (AResponses(i,e)==0)
          count(1)= count(1)+1;
      else
          count(2)= count(2)+1;
      end
   end
   
   %Determine the majority responses
   majority = 1;
   if (count(2)>count(1))
     majority =2;
   end
   
   %Determine whether the challenge has reliable 
   if (count(majority)<Threshold)
      Rate = Rate +1;
   end
     
end

end

