function [strf,options]=trnSCG(strf,datIdx,options,varargin)

%function [strf,options]=trnSCG(strf,datIdx,options)
%
% Scaled conjugate gradient optimization of STRF
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
% options = trnSCG
%   returns default options
%
%
% SEE ALSO: trnGradDesc, trnDirectFit
%
%(Some code modified from NETLAB)
% ====================


% Set default option values
% --------------------
global globDat;
optDeflt.funcName='trnSCG';
optRange.funcName={'trnSCG'};
optDeflt.display = 0;
optRange.display = [-Inf Inf];
optDeflt.maxIter = 5000;
optRange.maxIter = [0 Inf];
optDeflt.minErrChange = double(5*eps('single'));
optRange.minErrChange = double([eps('single') Inf]);
optDeflt.minStepSize = double(5*eps('single'));
optRange.minStepSize = double([eps('single') Inf]);

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
  set(gcf,'Units','points','Position',[20 50 800 400]);

    subplot(1,2,1);
    xlabel('parameter');
    ylabel('magnitude');
    subplot(1,2,2);
    xlabel('iteration');
    ylabel('estimation error');

end

% init
options.diagnostics.errs = [];  % record of estimation set error
h1 = []; h2 = []; h3 = [];      % plot handles

[strf, x] = strfPak(strf);              % initial seed

estbadcnt = 0;
stopbadcnt = 0;


nparams = length(x);

sigma0 = 1.0e-4;
[strf, fold] = strfErr(strf,datIdx);	% Initial function value.
fnow = fold;

options.diagnostics.errs = fnow;

[strf, gradnew] = strfGrad(strf,datIdx); % Initial gradient.
gradold = gradnew;

d = -gradnew;				% Initial search direction.
success = 1;				% Force calculation of directional derivs.
nsuccess = 0;				% nsuccess counts number of successes.
beta = 1.0;				% Initial scale parameter.
betamin = 1.0e-15; 			% Lower bound on scale.
betamax = 1.0e100;			% Upper bound on scale.
j = 1;					% j counts number of iterations.

etime = 0;
% Main optimization loop.
while (j <= options.maxIter)

    tic;
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
      return
    end
    sigma = sigma0/sqrt(kappa);
    xplus = x + sigma*d;
    [strf, gplus] = strfGrad(strfUnpak(strf,xplus),datIdx);
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

  options.diagnostics.errs=[options.diagnostics.errs fnow];

  if options.display < 0
    if j==1 || mod(j,options.display)==0
    
        delete([h1 h2]);
        subplot(1,2,1); hold on; %axis tight;
        h1 = bar(x(1:end-1)); title('Parameter magnitudes');
        subplot(1,2,2); hold on;
        h2 = plot(options.diagnostics.errs); title('Training Set Error');
        drawnow;
    end
  elseif options.display > 0
    if j==1 || mod(j,options.display)==0
      fprintf('iter: %03d | error: %.3f | step: %0.4f | etime: %0.0fs\n',j,fnow,max(abs(alpha*d)), etime);
    end
  end

  if (success == 1)
    % Test for termination

    if (max(abs(alpha*d)) < options.minStepSize && max(abs(fnew-fold)) < options.minErrChange)
      fprintf('minErrChange reached.\n');
      strf=strfUnpak(strf,x);
      return;

    else
      % Update variables for new position
      fold = fnew;
      gradold = gradnew;
      [strf,gradnew] = strfGrad(strfUnpak(strf,x),datIdx);;
     
      % If the gradient is zero then we are done.
      if (gradnew*gradnew' == 0)
	strf=strfUnpak(strf,x);
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
  etime = toc;
end

% If we get here, then we haven't terminated in the given number of 
% iterations.

strf=strfUnpak(strf,x);
