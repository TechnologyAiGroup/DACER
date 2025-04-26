function [AResponse] = ComputeResponseMXPUF2( ...
                       x_XPw,y_XPw,x,y,fb_a, ...
                       AChallenge,nChallenge,chalSize ...
                       )
      
      AResponse = zeros(nChallenge,1);
                      
      %Compute the Phi_xXORPUF based on AChallenge,nChallengs,chalSize and Transform
      %function
      chalX = chalSize;
      Phi_xXORPUF = Transform(AChallenge,nChallenge,chalX);
      %Compute the response of xXORPUF based on Phi_xXORPUF,xnewXORw and
      %ComputeResponseXOR function with Size = chalX+1
      Phi_xXORPUF = fliplr(Phi_xXORPUF);
      [response_xXORPUF, ~] = classify(Phi_xXORPUF,x_XPw,x);
             
      chalY = chalSize+1;
      %Create the challenge set for yXOR PUF
      yXORPUF_MXPUF_MSet = zeros(nChallenge,chalY);
      for r=1:nChallenge
          for c=1:(fb_a-1)
              yXORPUF_MXPUF_MSet(r,c) = AChallenge(r,c);
          end
      end 

      for r=1:nChallenge
          for c=(fb_a+1):chalY
              yXORPUF_MXPUF_MSet(r,c) = AChallenge(r,c-1);
          end
      end   
      
      for r=1:nChallenge
          yXORPUF_MXPUF_MSet(r,fb_a)=response_xXORPUF(r);
      end
      %Compute the Phi_yXORPUF based on yXORPUF_MSet,nMSet,chalY and Transform
      %function
      Phi_yXORPUF_MXPUF = Transform(yXORPUF_MXPUF_MSet,nChallenge,chalY);
      %Compute the response of yXORPUF based on Phi_yXORPUF,ynewXORw and
      %ComputeResponseXOR function with Size = chalY+1
      response_yXORPUF_MXPUF = ComputeResponseXOR(y_XPw,y, ... 
                               Phi_yXORPUF_MXPUF,nChallenge,chalY+1);      
      
      %In case of MXPUF, the response of y-XOR PUF is the output of MXPUF
      for r=1:nChallenge
         AResponse(r) = response_yXORPUF_MXPUF(r);
      end      
                   
end

