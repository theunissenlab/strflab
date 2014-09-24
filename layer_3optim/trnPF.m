function [strf,options]=trnPF(strf,datIdx,options,varargin)

%function [strf,options]=trnPF(strf,datIdx,options)
%
% Pathfinding optimization of STRF with arbitrary gaussian prior, determines best weight for
% prior using stopping set
%
% INPUT:
%            [strf] = model structure obtained via upper level *Init functions
%          [datIdx] = a vector containing indices of the samples to be used in the 
%                     fitting.
%         [options] = option structure containing fields:
%         .funcName = char array, name of this function.
%          .display = Number of iterations between displaying of errors as we train.
%                     Positive integer means display in stdout.
%                     Negative integer means display in a figure window.
%                     0 means do not display errors (default)
%          .maxIter = Maximum number of training iterations (default: 5000)
%     .minErrChange = The (positive real) size of error change which determines convergence
%	                  (default = eps(class(stim)))
%      .minStepSize = The (positive real) size of step size which determines convergence
%	                  (default = eps(class(stim)))
%        .lamdaInit = Initial weight on prior, will be decreased after each iteration of fitting (default: 2000000)
%         .lamdaDiv = Value by which current lambda is divided after each iteration of fitting (default:2)
%                .A = Covariance matrix for Gaussian prior (i.e. identity matrix for ridge regression)
%
%
% OUTPUT:
%            [strf] = structure containing model fit by trnGradDesc.
%         [options] = option structure with additional fields:
%      .diagnostics = contains information on the fitting process.
% .diagnostics.errs = this has dimensions 1 x N where this has the error on the
%                     training set, over the course of the fitting iterations.
%
% EXAMPLE:
% options = trnPF
%   returns default options
%
%
% SEE ALSO: trnGradDesc, trnSCG
%
%(Some code modified from NETLAB)
% ====================


% Set default option values
% --------------------
global globDat;
optDeflt.funcName='trnPF';
optRange.funcName={'trnPF'};
optDeflt.display = 0;
optRange.display = [-Inf Inf];
optDeflt.maxIter = 5000;
optRange.maxIter = [0 Inf];
optDeflt.minErrChange = double(5*eps('single'));
optRange.minErrChange = double([eps('single') Inf]);
optDeflt.minStepSize = double(5*eps('single'));
optRange.minStepSize = double([eps('single') Inf]);
optDeflt.is1D = 0;
optRange.is1D = [0 1];
optDeflt.nDispTruncate = 2;
optRange.nDispTruncate = [0 Inf];
optDeflt.errLastN=50;
optRange.errLastN=[-Inf Inf];
optDeflt.lamdaInit=2000000;
optRange.lamdaInit=[1 Inf];
optDeflt.lamdaDiv=2;
optRange.lamdaDiv=[1.01 Inf];
if nargin<3
  options=optDeflt;
else
  options=defaultOpt(options,optDeflt,optRange);
end
if nargin<1
  strf=optDeflt;
  return;
end


% initialize plot
if options.display < 0
  clf;
  set(gcf,'Units','points','Position',[20 50 1000 300]);
  subplot(1,3,1);
  if options.is1D == 0
      xlabel('parameter');
      ylabel('magnitude');
  end
  subplot(1,3,2);
  xlabel('iteration');
  ylabel('estimation error');
  subplot(1,3,3);
  xlabel('iteration');
  ylabel('stopping error');

end

% init
options.diagnostics.errs = [];  % record of estimation set error
h1 = []; h2 = []; h3 = [];      % plot handles

[strf, x] = strfPak(strf);              % initial seed

estbadcnt = 0;
stopbadcnt = 0;


nparams = length(x);

%% Set A matrix
A=options.A;

sigma0 = 1.0e-4;
[strf, fold] = strfErr(strf,datIdx);  % Initial function value.
[strf, besterr] = strfErr(strf,varargin{:});

% iter = 1;
j = 1;                    % j counts number of iterations.
regerr=x*(A*x')/2;
lam=options.lamdaInit*fold/(regerr+1);

manylog=[];

while lam>.1

manylog=[manylog;x];

lam=lam/options.lamdaDiv;

[strf, fold] = strfErr(strf,datIdx);	% Initial function value.
%%
regerr=x*(A*x')/2;
fold=fold+regerr*lam;

fnow = fold;

% options.diagnostics.errs = fnow;

[strf, gradnew] = strfGrad(strf,datIdx); % Initial gradient.
gradnew = gradnew+lam*x*A;
gradold = gradnew;

d = -gradnew;				% Initial search direction.
success = 1;				% Force calculation of directional derivs.
nsuccess = 0;				% nsuccess counts number of successes.
beta = 1.0;				% Initial scale parameter.
betamin = 1.0e-15; 			% Lower bound on scale.
betamax = 1.0e100;			% Upper bound on scale.
% j = 1;                    % j counts number of iterations.
flog(j, :) = fold;
pointlog(j, :) = x;



% Main optimization loop.
while (j <= options.maxIter)

  % Calculate first and second directional derivatives.
  if (success == 1)
    mu = d*gradnew';
    if (mu >= 0)
      d = - gradnew;
      mu = d*gradnew';
    end
    kappa = d*d';
    if kappa < eps
      % Ryan what is this please fix...
      %options(8) = fnow;
      keyboard
      return
    end
    sigma = sigma0/sqrt(kappa);
    xplus = x + sigma*d;
    [strf, gplus] = strfGrad(strfUnpak(strf,xplus),datIdx);
    gplus = gplus+lam*xplus*A;
    theta = (d*(gplus' - gradnew'))/sigma;
  end

  % Increase effective curvature and evaluate step size alpha.
  delta = theta + beta*kappa;
  if (delta <= 0) 
    delta = beta*kappa;
    beta = beta - theta/kappa;
  end
  alpha = - mu/delta;
  
  % Calculate the comparison ratio.
  xnew = x + alpha*d;
  [strf, fnew]= strfErr(strfUnpak(strf,xnew),datIdx);
  fnew = fnew+lam*xnew*(A*xnew')/2;
  
  Delta = 2*(fnew - fold)/(alpha*mu);
  if (Delta  >= 0)
    success = 1;
    nsuccess = nsuccess + 1;
    x = xnew;
    fnow = fnew;
  else
    success = 0;
    fnow = fold;
  end


    % Store relevant variables
    flog(j) = fnow;		% Current function value

    pointlog(j,:) = x;	% Current position

	scalelog(j) = beta;	% Current scale parameter


  % options.diagnostics.errs=[options.diagnostics.errs fnow];
  options.diagnostics.errs(1,j) = fnow;
  [strf,stoperr] = strfErr(strf,varargin{:});
  options.diagnostics.errs(2,j) = stoperr;

  if options.display < 0
    if j==1 || mod(j,options.display)==0
    
        delete([h1 h2 h3]);
        subplot(1,3,1); hold on; %axis tight;
        if options.is1D
            simagesc(reshape(x(1:(end-1)),strf.nIn,length(strf.delays)));
            axis tight;
        else
            h1 = bar(x(1:end-1));
        end
        
        subplot(1,3,2); hold on;
        h2 = plot((1+optDeflt.nDispTruncate):length(options.diagnostics.errs),options.diagnostics.errs(1,(1+optDeflt.nDispTruncate):end),'-');
        subplot(1,3,3); hold on;
        h3 = plot((1+optDeflt.nDispTruncate):length(options.diagnostics.errs),options.diagnostics.errs(2,(1+optDeflt.nDispTruncate):end),'-');
        drawnow;
    end
  elseif options.display > 0
    if j==1 || mod(j,options.display)==0
      fprintf('iter: %03d | error: %.3f Scale %e\n',j,fnow,beta);
    end
  end

  if (success == 1)
    % Test for termination

    if stoperr < besterr
      stopbadcnt = 0;
      besterr = stoperr;
      bestx = x;
    else
      stopbadcnt = stopbadcnt + 1;
    end

    if (options.errLastN > 0 && stopbadcnt >= options.errLastN)
        % j = j + 1;
        % iter = iter + 1;
        % stopbadcnt = 0;
        % break;
        strf=strfUnpak(strf,bestx);
        strf.flog = flog;
        strf.pointlog = pointlog;
        strf.scalelog = scalelog;
        break
    end

    if (max(abs(alpha*d)) < options.minStepSize & max(abs(fnew-fold)) < options.minErrChange)
      % strf=strfUnpak(strf,x);
      % [strf,stoperr] = strfErr(strf,varargin{:});
      % if stoperr < besterr
      %     besterr = stoperr;
      %     bestx = x;
      % end
      % options.diagnostics.errs(2,iter) = stoperr;
      j = j + 1;
      % iter = iter + 1;
      % stopbadcnt = 0;
      break;

    else
      % Update variables for new position
      fold = fnew;
      gradold = gradnew;
      [strf,gradnew] = strfGrad(strfUnpak(strf,x),datIdx);
      gradnew = gradnew + lam*x*A;
     
      % If the gradient is zero then we are done.
      if (gradnew*gradnew' == 0)
	      strf=strfUnpak(strf,x);
	    keyboard
	    return;
      end
    end
  end

  % Adjust beta according to comparison ratio.
  if (Delta < 0.25)
    beta = min(4.0*beta, betamax);
  end
  if (Delta > 0.75)
    beta = max(0.5*beta, betamin);
  end

  % Update search direction using Polak-Ribiere formula, or re-start 
  % in direction of negative gradient after nparams steps.
  if (nsuccess == nparams)
    d = -gradnew;
    nsuccess = 0;
  else
    if (success == 1)
      gamma = (gradold - gradnew)*gradnew'/(mu);
      d = gamma*d - gradnew;
    end
  end
  j = j + 1;
  % iter = iter + 1;
end

end

% If we get here, then we haven't terminated in the given number of 
% iterations.

strf=strfUnpak(strf,bestx);
strf.flog = flog;
strf.pointlog = pointlog;
strf.scalelog = scalelog;