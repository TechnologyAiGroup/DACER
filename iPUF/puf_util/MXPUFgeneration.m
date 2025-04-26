function [x_XPw,y_XPw] = MXPUFgeneration(x,y,chalSize,mu,sigma)

  
  x_XPw = normrnd(mu,sigma,x,chalSize+1);
  y_XPw = normrnd(mu,sigma,y,chalSize+2);
  
end

