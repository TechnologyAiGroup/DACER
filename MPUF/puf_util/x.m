function [xopt,fopt] = x(f,d,x)

% Written by Sun Yi (yi@idsia.ch).

% parameters
L = 4+3*floor(log(d));
etax = 1; etaA = 0.5*min(1.0/d,0.25);
shape = max(0.0, log(L/2+1.0)-log(1:L)); shape = shape / sum(shape);

% initialize
xopt = x; 
%fopt = f(x);
fopt= feval(f, x);
A = zeros(d);
weights = zeros(1,L);
fit = zeros(1,L);


while 0<1
    expA = expm(A);
    
    % step 1: sampling & importance mixing
    Z = randn(d,L); X = repmat(x,1,L)+expA*Z;
    for i = 1 : L, fit(i) = f(X(:,i)); end

    % step 2: fitness reshaping
    [~, idx] = sort(fit); weights(idx) = shape;
    if fit(idx(1)) < fopt
        xopt = X(:,idx(1)); fopt = fit(idx(1));
    end

    % step 3: compute the gradient for C and x
    G = (repmat(weights,d,1).*Z)*Z' - sum(weights)*eye(d);
    dx = etax * expA * (Z*weights');
    dA = etaA * G;
  
    % step 4: compute the update  
    x = x + dx; A = A + dA;

    if trace(A)/d < -10*log(10), break; end
end