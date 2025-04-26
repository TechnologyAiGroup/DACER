function [AResponse] = ComputeNoisyResponsesMXPUF(x_XPw,y_XPw,x,y,fb_a,...
                                     AChallenge,nChallenge,chalSize,mu,sigma,sigmaNoise,Evaluations)

 AResponse = zeros(Evaluations,nChallenge);

 for e=1:Evaluations  
     
      %Create a noise vector 
      sigNoise = sigmaNoise*sigma;
      
      %generate (x,y)-MXPUF

      [x_XPw_Noise,y_XPw_Noise]=MXPUFgeneration(x,y,chalSize,mu,sigNoise);
      
      %Create a MXPUF affected by the noise [x_XPw_Noise,y_XPw_Noise]
      new_x_XPw = x_XPw + x_XPw_Noise;
      new_y_XPw = y_XPw + y_XPw_Noise;
      
      %Evaluate the responses
    
      AResponse(e,:)= ComputeResponseMXPUF( ...
                                           new_x_XPw,new_y_XPw,x,y,fb_a, ...
                                           AChallenge,nChallenge,chalSize ...
                                           );
            

      
  end  
end

