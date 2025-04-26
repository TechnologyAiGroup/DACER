function [ARank,AOccurence] = ComputeRankingModel(Matching, Noise, nRows, nCols)
% This function computes the rank of found model wrt noise rate.
% nRows and nCols means the number of runs of CMAES and the number of APUFs in XORPUF
%Matching[i,j]=1/0 means whether the CMAES can/ cannot see the model Aj in run i 
%Noise[i,j]  shows the noise rate of Aj in run i 

%We also have an array ACoccurence to measure the frequency


  
  ARank = zeros(nRows,1);
  %The rank should be from 1 to nXOR= nCols
  
  for i=1:nRows 
      Seeflag = 0; % Flag shows we can see any model in run i or not
      position   = -1;% If we can see, then what is the APUF 
      
      for j=1:nCols %A1, ..,AnXOR
          if(Matching(i,j)==1)
              Seeflag = 1;
              position = j;
              break;
          end              
      end
      

      
      if(Seeflag==1)
         ARank(i) = 1;
         Value = Noise(i,position);
         for j=1:nCols
             if((Noise(i,j)>Value)&&(j~=position))
                 ARank(i) = ARank(i)+1;
             end
         end
      end
      
  end 
  
  AOccurence = zeros(nCols+1,1);
  
  for i=1:nRows
      AOccurence(ARank(i)+1) =  AOccurence(ARank(i)+1)+1;
  end

  for i=1:(nCols+1)
      AOccurence(i) =  AOccurence(i)/nRows;
  end


%   er  = 0 ; 
%   m = size(AChallenge,1);
%   for i=1:m
%       R = apuf.getResponse(AChallenge(i,:));
%       if(R ~= AResponse(i))
%          er = er+1;
%       end
%   end
% 
%    error = er/m; 
end

