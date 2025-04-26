function [xmin, ...      % minimum search point of last iteration
	  fmin, ...      % function value of xmin
	  counteval, ... % number of function evaluations done
	  stopflag, ...  % stop criterion reached
	  out, ...     % struct with various histories and solutions
	  bestever ... % struct containing overall best solution (for convenience)
	 ] = cmaes( ...
    fitfun, ...    % name of objective/fitness function
    xstart, ...    % objective variables initial point, determines N
    insigma, ...   % initial coordinate wise standard deviation(s)
    inopts, ...    % options struct, see defopts below
    x, ...
    h, ...
    flag_diag,flag_earlystop,...
    flgFitnessLoss, loss_alpha, loss_func, dis_mu, dis_sig, ...
    AResponse_TrS1,nXOR1,Size1,AResponse_TrS_final1,randa1,Phi_TeS1_1,...
    varargin)     % arguments passed to objective function 

definput.fitfun = 'modelAcc'; % frosen; fcigar; see end of file for more
definput.xstart = rand(10,1); % 0.50*ones(10,1);
definput.sigma = 0.5;
global x1;
global h1;
global nXOR;
global Size;
global randa;
global AResponse_TrS;
global AResponse_TrS_final;
global Phi_TrS1;
x1=x;
h1=h;
nXOR = nXOR1;
Size = Size1;
randa = randa1;
AResponse_TrS = AResponse_TrS1;
AResponse_TrS_final = AResponse_TrS_final1;
Phi_TrS1=Phi_TeS1_1;
flag=0;
pr_re=0;

defopts.StopFitness  = '1e-10 % stop if f(xmin) < stopfitness, minimization';
defopts.MaxFunEvals  = '100000  % maximal number of fevals';
defopts.MaxIter      = '1e3*(N+5)^2/sqrt(popsize) % maximal number of iterations';
defopts.Restarts     = '0    % number of restarts ';
defopts.IncPopSize   = '2    % multiplier for population size before each restart';
defopts.PopSize      = '(4 + floor(3*log(N)))  % population size, lambda'; 
defopts.ParentNumber = 'floor(popsize/2)       % AKA mu, popsize equals lambda';
defopts.RecombinationWeights = 'superlinear decrease % or linear, or equal';
defopts.DiagonalOnly = '0  % C is diagonal for given iterations, 1==always';
defopts.CMA.cs = '(mueff+2)/(N+mueff+3)  % cumulation constant for step-size'; 
defopts.CMA.damps = '1 + 2*max(0,sqrt((mueff-1)/(N+1))-1) + cs  % damping for step-size';
defopts.CMA.ccum = '(4 + mueff/N) / (N+4 + 2*mueff/N)  % cumulation constant for pc';
defopts.CMA.ccov1 = '2 / ((N+1.3)^2+mueff)  % learning rate for rank-one update'; 
defopts.CMA.ccovmu = '2 * (mueff-2+1/mueff) / ((N+2)^2+mueff) % learning rate for rank-mu update'; 
defopts.CMA.on     = 'yes'; 
defopts.CMA.active = '0  % active CMA 1: neg. updates with pos. def. check, 2: neg. updates'; 
defopts.DiagonalOnly=flag_diag;
  
% ---------------------- Handling Input Parameters ----------------------

input.fitfun = fitfun; % record used input
input.xstart = xstart;
input.sigma = insigma;
sigma=insigma;
% Compose options opts
opts = defopts;
counteval = 0; countevalNaN = 0; 
irun = 0;
while irun <= myeval(opts.Restarts) % for-loop does not work with resume
  irun = irun + 1; 

% ------------------------ Initialization -------------------------------

% Handle resuming of old run
xmean = myeval(xstart); 
if all(size(xmean) > 1)
   xmean = mean(xmean, 2); % in case if xstart is a population
elseif size(xmean, 2) > 1
  xmean = xmean';
end 

  % Assign settings from input parameters and options for myeval...
  N = size(xmean, 1);  
  lambda0 = floor(myeval(opts.PopSize) * myeval(opts.IncPopSize)^(irun-1)); 
  % lambda0 = floor(myeval(opts.PopSize) * 3^floor((irun-1)/2)); 
  popsize = lambda0;
  lambda = lambda0;
  insigma = myeval(insigma);
  if all(size(insigma) == [N 2]) 
    insigma = 0.5 * (insigma(:,2) - insigma(:,1));
  end

  
%--------------------------------------------------------------
% Evaluate options
flgDiagonalOnly = myeval(opts.DiagonalOnly); 
MaxFunEvals = myeval(opts.MaxFunEvals);
flgActiveCMA = myeval(opts.CMA.active); 
bnd.isactive=0;
  %countiter = countiter + 1; 
chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));
  fitness.hist=NaN*ones(1,10+ceil(3*10*N/lambda)); % history of fitness values
  fitness.histsel=NaN*ones(1,10+ceil(3*10*N/lambda)); % history of fitness values
  fitness.histbest=[]; % history of fitness values
  fitness.histmedian=[]; % history of fitness values
    fitness.hist(1)=feval(fitfun, xmean, varargin{:}); 
    %fitness.hist(1)=fitness.hist(1)-0.1*feval("model_corr", xmean, varargin{:});%new
    fitness.histsel(1)=fitness.hist(1);


% -------------------- Generation Loop --------------------------------

countiter = 0;
  sigma = max(insigma);              % overall standard deviation
  pc = zeros(N,1); ps = zeros(N,1);  % evolution paths for C and sigma

  if length(insigma) == 1
    insigma = insigma * ones(N,1) ;
  end
  diagD = insigma/max(insigma);      % diagonal matrix D defines the scaling
  diagC = diagD.^2; 
  if flgDiagonalOnly ~= 1            % use at some point full covariance matrix
    B = eye(N,N);                      % B defines the coordinate system
    BD = B.*repmat(diagD',N,1);        % B*D for speed up only
    C = diag(diagC);                   % covariance matrix == BD*(BD)'
  end
  if flgDiagonalOnly
    B = 1; 
  end
while (counteval<MaxFunEvals)
  % set internal parameters
  if countiter == 0
    % Strategy internal parameter setting: Selection  
    mu = myeval(opts.ParentNumber); % number of parents/points for recombination
    if strncmp(lower(opts.RecombinationWeights), 'equal', 3)
      weights = ones(mu,1); % (mu_I,lambda)-CMA-ES
    elseif strncmp(lower(opts.RecombinationWeights), 'linear', 3)
      weights = mu+0.5-(1:mu)'; 
    elseif strncmp(lower(opts.RecombinationWeights), 'superlinear', 3)
      % use (lambda+1)/2 as reference if mu < lambda/2
      weights = log(max(mu, lambda/2) + 1/2)-log(1:mu)'; % muXone array for weighted recombination
    else
      error(['Recombination weights to be "' opts.RecombinationWeights ...
             '" is not implemented']);
    end
    mueff=sum(weights)^2/sum(weights.^2); % variance-effective size of mu
    weights = weights/sum(weights);     % normalize recombination weights array
    if mueff == lambda
      error(['Combination of values for PopSize, ParentNumber and ' ...
             ' and RecombinationWeights is not reasonable']);
    end
    
    % Strategy internal parameter setting: Adaptation
    cc = myeval(opts.CMA.ccum); % time constant for cumulation for covariance matrix
    cs = myeval(opts.CMA.cs); 

    if myevalbool(opts.CMA.on) 
      ccov1 = myeval(opts.CMA.ccov1); 
      ccovmu = min(1-ccov1, myeval(opts.CMA.ccovmu));
    else
      ccov1 = 0;
      ccovmu = 0;
    end
    
    % flgDiagonalOnly = -lambda*4*1/ccov; % for ccov==1 it is not needed
    % 0 : C will never be diagonal anymore
    % 1 : C will always be diagonal
    % >1: C is diagonal for first iterations, set to 0 afterwards
    if flgDiagonalOnly < 1
      flgDiagonalOnly = 0; 
    end
    if flgDiagonalOnly
      ccov1_sep = min(1, ccov1 * (N+1.5) / 3); 
      ccovmu_sep = min(1-ccov1_sep, ccovmu * (N+1.5) / 3); 
    end
 
    damps = myeval(opts.CMA.damps); 

  end    

  % Generate and evaluate lambda offspring
 
  fitness.raw = repmat(NaN, 1, lambda);


  for k=find(isnan(fitness.raw)) 
    % fitness.raw(k) = NaN; 

    % Resample, until fitness is not NaN
    while isnan(fitness.raw(k))
      if k <= lambda  % regular samples (not the re-evaluation-samples)
        arz(:,k) = randn(N,1); % (re)sample

        if flgDiagonalOnly  
          arx(:,k) = xmean + sigma * diagD .* arz(:,k);              % Eq. (1)
        else
          arx(:,k) = xmean + sigma * (BD * arz(:,k));                % Eq. (1)
        end
      else % re-evaluation solution with index > lambda
        if flgDiagonalOnly  
          arx(:,k) = arx(:,k-lambda) + (noiseEpsilon * sigma) * diagD .* randn(N,1);
        else
          arx(:,k) = arx(:,k-lambda) + (noiseEpsilon * sigma) * (BD * randn(N,1));
        end
      end
      
        arxvalid(:,k) = arx(:,k);

      fitness.raw(k) = feval(fitfun, arxvalid(:,k), varargin{:}); 
      fitness.raw(k) = fitness.raw(k)- 0.1*feval("model_corr", arxvalid(:,k), varargin{:});
      fitness.raw(k) = fitness.raw(k)- 30*feval("model_bin", arxvalid(:,k), varargin{:});
      if flgFitnessLoss
        fit_loss = feval(loss_func, fitness.raw(k), dis_mu, dis_sig);          
        fitness.raw(k) = fitness.raw(k) - loss_alpha * fit_loss;
      end
      if isnan(fitness.raw(k))
	countevalNaN = countevalNaN + 1;
      end
    end
    counteval = counteval + 1; % retries due to NaN are not counted
  end

 
  fitness.sel = fitness.raw; 
  
  % Sort by fitness 
  [fitness.raw, fitness.idx] = sort(fitness.raw,'descend'); 
  [fitness.sel, fitness.idxsel] = sort(fitness.sel,'descend');  % minimization
  fitness.hist(2:end) = fitness.hist(1:end-1);    % record short history of
  fitness.hist(1) = fitness.raw(1);               % best fitness values
  if length(fitness.histbest) < 120+ceil(30*N/lambda) || ...
       (mod(countiter, 5) == 0  && length(fitness.histbest) < 2e4)  % 20 percent of 1e5 gen.
    fitness.histbest = [fitness.raw(1) fitness.histbest];          % best fitness values
    fitness.histmedian = [median(fitness.raw) fitness.histmedian]; % median fitness values
  else
    fitness.histbest(2:end) = fitness.histbest(1:end-1); 
    fitness.histmedian(2:end) = fitness.histmedian(1:end-1); 
    fitness.histbest(1) = fitness.raw(1);           % best fitness values
    fitness.histmedian(1) = median(fitness.raw);    % median fitness values
  end
  fitness.histsel(2:end) = fitness.histsel(1:end-1); % record short history of
  fitness.histsel(1) = fitness.sel(1);               % best sel fitness values

  if flag_earlystop==1
    if fitness.sel(1)-pr_re<10^(-6)
        flag=flag+1;
    else
        flag=0;
    end
    if flag>=10
        break;
    end
    pr_re=fitness.sel(1);
  end
  
  % Calculate new xmean, this is selection and recombination 
  xold = xmean; % for speed up of Eq. (2) and (3)
  cmean = 1;  % 1/min(max((lambda-1*N)/2, 1), N);  % == 1/kappa
  xmean = (1-cmean) * xold + cmean * arx(:,fitness.idxsel(1:mu))*weights; 
  zmean = arz(:,fitness.idxsel(1:mu))*weights;%==D^-1*B'*(xmean-xold)/sigma
  if mu == 1
    fmean = fitness.sel(1);
  else
    fmean = NaN; % [] does not work in the latter assignment
  end
  
  % Cumulation: update evolution paths
  ps = (1-cs)*ps + sqrt(cs*(2-cs)*mueff) * (B*zmean);          % Eq. (4)
  hsig = norm(ps)/sqrt(1-(1-cs)^(2*countiter))/chiN < 1.4 + 2/(N+1);

  pc = (1-cc)*pc ...
        + hsig*(sqrt(cc*(2-cc)*mueff)/sigma/cmean) * (xmean-xold);     % Eq. (2)


  % Adapt covariance matrix
  neg.ccov = 0;  % TODO: move parameter setting upwards at some point
  if ccov1 + ccovmu > 0                                                    % Eq. (3)
    if flgDiagonalOnly % internal linear(?) complexity
      diagC = (1-ccov1_sep-ccovmu_sep+(1-hsig)*ccov1_sep*cc*(2-cc)) * diagC ... % regard old matrix 
          + ccov1_sep * pc.^2 ...               % plus rank one update
          + ccovmu_sep ...                      % plus rank mu update
            * (diagC .* (arz(:,fitness.idxsel(1:mu)).^2 * weights));
%             * (repmat(diagC,1,mu) .* arz(:,fitness.idxsel(1:mu)).^2 * weights);
      diagD = sqrt(diagC); % replaces eig(C)
    else
      arpos = (arx(:,fitness.idxsel(1:mu))-repmat(xold,1,mu)) / sigma;
      % "active" CMA update: negative update, in case controlling pos. definiteness 
      if flgActiveCMA > 0
        % set parameters
        neg.mu = mu;  
        neg.mueff = mueff;
        if flgActiveCMA > 10  % flat weights with mu=lambda/2
          neg.mu = floor(lambda/2);  
          neg.mueff = neg.mu;
        end

        neg.ccov = (1 - ccovmu) * 0.25 * neg.mueff / ((N+2)^1.5 + 2*neg.mueff);
        neg.minresidualvariance = 0.66;  % keep at least 0.66 in all directions, small popsize are most critical
        neg.alphaold = 0.5;  % where to make up for the variance loss, 0.5 means no idea what to do
                             % 1 is slightly more robust and gives a better "guaranty" for pos. def., 
                             % but does it make sense from the learning perspective for large ccovmu? 

        neg.ccovfinal = neg.ccov;

        % prepare vectors, compute negative updating matrix Cneg and checking matrix Ccheck
        arzneg = arz(:,fitness.idxsel(lambda:-1:lambda - neg.mu + 1));
        [arnorms idxnorms] = sort(sqrt(sum(arzneg.^2, 1))); 
        [ignore idxnorms] = sort(idxnorms);  % inverse index 
        arnormfacs = arnorms(end:-1:1) ./ arnorms; 
        arnorms = arnorms(end:-1:1); % for the record
        if flgActiveCMA < 20
          arzneg = arzneg .* repmat(arnormfacs(idxnorms), N, 1);  % E x*x' is N
        end
        if flgActiveCMA < 10 && neg.mu == mu  % weighted sum
          if mod(flgActiveCMA, 10) == 1 % TODO: prevent this with a less tight but more efficient check (see below) 
            Ccheck = arzneg * diag(weights) * arzneg';  % in order to check the largest EV
          end
          artmp = BD * arzneg;
          Cneg = artmp * diag(weights) * artmp';
        else  % simple sum
          if mod(flgActiveCMA, 10) == 1
            Ccheck = (1/neg.mu) * arzneg*arzneg';  % in order to check largest EV
          end
          artmp = BD * arzneg;
          Cneg = (1/neg.mu) * artmp*artmp';

        end

        % check pos.def. and set learning rate neg.ccov accordingly, 
        % this check makes the original choice of neg.ccov extremly failsafe 
        % still assuming C == BD*BD', which is only approxim. correct 
        if mod(flgActiveCMA, 10) == 1 && 1 - neg.ccov * arnorms(idxnorms).^2 * weights < neg.minresidualvariance
          maxeigenval = max(eig(Ccheck));  % norm is much slower, because (norm()==max(svd())
          neg.ccovfinal = min(neg.ccov, (1-ccovmu)*(1-neg.minresidualvariance)/maxeigenval); 
                                        % -ccov1 removed to avoid error message??
        end
        % update C
        C = (1-ccov1-ccovmu+neg.alphaold*neg.ccovfinal+(1-hsig)*ccov1*cc*(2-cc)) * C ... % regard old matrix 
            + ccov1 * pc*pc' ...     % plus rank one update
            + (ccovmu + (1-neg.alphaold)*neg.ccovfinal) ...  % plus rank mu update
              * arpos * (repmat(weights,1,N) .* arpos') ...
              - neg.ccovfinal * Cneg;                        % minus rank mu update
      else  % no active (negative) update
        C = (1-ccov1-ccovmu+(1-hsig)*ccov1*cc*(2-cc)) * C ... % regard old matrix 
            + ccov1 * pc*pc' ...     % plus rank one update
            + ccovmu ...             % plus rank mu update
              * arpos * (repmat(weights,1,N) .* arpos');
      end
      diagC = diag(C);
    end
  end
  
  % Adapt sigma

    % exp(1) is still not reasonably small enough
    sigma = sigma * exp(min(1, (sqrt(sum(ps.^2))/chiN - 1) * cs/damps));   

  % Update B and D from C

  if ~flgDiagonalOnly && (ccov1+ccovmu+neg.ccov) > 0 && mod(countiter, 1/(ccov1+ccovmu+neg.ccov)/N/10) < 1
    C=triu(C)+triu(C,1)'; % enforce symmetry to prevent complex numbers
    [B,tmp] = eig(C);     % eigen decomposition, B==normalized eigenvectors
                          % effort: approx. 15*N matrix-vector multiplications
    diagD = diag(tmp); 

    diagC = diag(C); 
    diagD = sqrt(diagD); % D contains standard deviations now
    BD = B.*repmat(diagD',N,1); % O(n^2)
  end % if mod

  % Align/rescale order of magnitude of scales of sigma and C for nicer output
  if 1 < 2 && sigma > 1e10*max(diagD) && sigma > 8e14 * max(insigma)
    fac = sigma; % / max(diagD);
    sigma = sigma/fac;
    pc = fac * pc;
    diagD = fac * diagD; 
    if ~flgDiagonalOnly
      C = fac^2 * C; % disp(fac);
      BD = B .* repmat(diagD',N,1); % O(n^2), but repmat might be inefficient todo?
    end
    diagC = fac^2 * diagC; 
  end

  if flgDiagonalOnly > 1 && countiter > flgDiagonalOnly 
    % full covariance matrix from now on 
    flgDiagonalOnly = 0; 
    B = eye(N,N);
    BD = diag(diagD);
    C = diag(diagC); % is better, because correlations are spurious anyway
  end


  % ----- numerical error management -----
  % Adjust maximal coordinate axis deviations

  % Adjust step size in case of (numerical) precision problem 
  if flgDiagonalOnly
    tmp = 0.1*sigma*diagD; 
  else
    tmp = 0.1*sigma*BD(:,1+floor(mod(countiter,N)));
  end
    
  
  %disp([num2str(counteval) ': ' num2str(fitness.hist(1))]);

  % ----- end output generation -----

end % while, end generation loop

% -------------------- Final Procedures -------------------------------

% Evaluate xmean and return best recent point in xmin
fmin = fitness.raw(1);
xmin = arxvalid(:, fitness.idx(1)); % Return best point of last generation.

end 

function res=myeval(s)
  if ischar(s)
    res = evalin('caller', s);
  else
    res = s;
  end
  
% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function res=myevalbool(s)
  if ~ischar(s) % s may not and cannot be empty
    res = s;
  else % evaluation string s
    if strncmpi(s, 'yes', 3) || strncmpi(s, 'on', 2) ...
	  || strncmpi(s, 'true', 4) || strncmp(s, '1 ', 2)
      res = 1;
    elseif strncmpi(s, 'no', 2) || strncmpi(s, 'off', 3) ...
	  || strncmpi(s, 'false', 5) || strncmp(s, '0 ', 2)
      res = 0;
    else
      try res = evalin('caller', s); catch
	error(['String value "' s '" cannot be evaluated']);
      end
      try res ~= 0; catch
	error(['String value "' s '" cannot be evaluated reasonably']);
      end
    end
  end
   
function res=myrange(x)
  res = max(x) - min(x);
  
 function f = modelAcc(w)
    global x1;
    global h1;
    global nXOR;
    global Size;
    global randa;
    global AResponse_TrS;
    r = 0;
    for i = 1:nXOR
        P = x1 * w((i-1)*(Size+1)+1:(i-1)*(Size+1) + Size,:);
        P = abs(P);
        for j=1:size(P,1)
        if P(j) > w((i-1)*(Size+1) + Size+1,:)
            P1(j,i)=1;
        else
            P1(j,i)=0;
        end
        end
    end
    for i=1:nXOR/2
       for j=1:size(P,1)
           if (P1(j,i)==1&&AResponse_TrS(i,j)==0)||(P1(j,randa(i))==1&&AResponse_TrS(randa(i),j)==0)||(P1(j,randa(i))==1&&AResponse_TrS(randa(i),j)==1&&P1(j,i)==1&&AResponse_TrS(i,j)==1)
               P_final(j,i)=1;
           else
               P_final(j,i)=0;
           end
       end
    end
    for i = 1:nXOR/2
        h_tilde = P_final(:,i);
        r = r + corr(h1',h_tilde,'Type','Pearson');
    end
    f= r;

function f = model_bin(w)
    global x1;
    global h1;
    global nXOR;
    global Size;
    global randa;
    global AResponse_TrS;
    global AResponse_TrS_final;
    global Phi_TrS1;
    r = 0;
    nTrS=size(x1,1);
    for i=1:nXOR
        w1(i,:)=w((i-1)*(Size+1)+1:(i-1)*(Size+1) + Size,:);
    end
    for i=1:nXOR
        AResponse(i,:) = ComputeResponseXOR(w1(i,:),1,x1,nTrS,Size);
    end
    for i=1:(nXOR/2)
        Response(i,:)=and(AResponse(i,:),AResponse(randa(i),:));
    end
    Response1(1,:)=Response(1,:);
    for i=2:(nXOR/2)
        Response1(1,:)=xor(Response1(1,:),Response(i,:));
    end
    r = crossentropy(Response1,AResponse_TrS_final(1,:));
    f= r;

 function f = model_corr(w)
    global nXOR;
    global Size;
    r = 0;
    for i = 1:nXOR-1
        for j = i+1:nXOR
        r = r + abs(corr(w((i-1)*Size+1:(i-1)*Size + Size,:),w((j-1)*Size+1:(j-1)*Size + Size,:),'Type','Pearson'));
        end
    end
    f= r;