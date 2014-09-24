function [strf, options] = trnSCG(strf,stim,resp,options,varargin)

%trnSCG	Scaled conjugate gradient optimization.
%
% [STRF, OPTIONS] = TRNSCG(STRF, STIM, RESP)
% [STRF, OPTIONS] = TRNSCG(STRF, STIM, RESP, OPTIONS)
% 	 OPTIONS  = TRNSCG
%
% Inputs:
%
% STRF (the seed)
%
% STIM
%
% RESP
%
% OPTIONS, A structure containing the training parameters
%
% The OPTIONS structure can include the following fields.  If a field
% is not included, the default value is used:
%

% OPTIONS.disp (integer)
%
%   Number of iterations between displaying of errors as we train.
%   If positive, display in stdout.
%   If negative, display in a figure window.
%   If 0, do not display errors (default, 0).
%
% OPTIONS.max_train_iter (nonnegative integer)
%   Maximum number of training iterations (default, 1000).
%   if .fixed_iter is supplied, this option has no effect.
%
%
% OPTIONS.err_change (positive real)
%	The size of error change which determines convergence
%	(default, eps)
%
% OPTIONS.min_step_size (positive real)
%	The size of step size which determines convergence
%	(default, eps)
%
%
% Outputs:
%
% STRF
%
% OPTIONS.diagnostics will contain information on the fitting process.
%
% OPTIONS.diagnostics.errs
%
%   This has dimensions 1 x N where this has the error on the
%   training set, over the course of the fitting iterations.
%	
%(Some code modified from NETLAB)

%  defaults

optDeflt.disp = 0;
optRange.disp = [-Inf Inf];
optDeflt.max_train_iter = 1000;
optRange.max_train_iter = [0 Inf];
optDeflt.err_change = eps;
optRange.err_change = [eps Inf];
optDeflt.min_step_size = eps;
optRange.min_step_size = [eps Inf];

if nargin < 4
  options = optDeflt;
else
  options = defaultOpt(options,optDeflt,optRange);
end
if nargin < 3
  strf = optDeflt;
  return;
end


% initialize plot
if options.disp < 0
  figure;
  set(gcf,'Units','points','Position',[100 100 800 400]);

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

x = strfPak(strf);              % initial seed

estbadcnt = 0;
stopbadcnt = 0;


nparams = length(x);

sigma0 = 1.0e-4;
fold = strfErr(strf,stim,resp);	% Initial function value.
fnow = fold;

options.diagnostics.errs = fnow;

gradnew = strfGrad(strf,stim,resp); % Initial gradient.
gradold = gradnew;

d = -gradnew;				% Initial search direction.
success = 1;				% Force calculation of directional derivs.
nsuccess = 0;				% nsuccess counts number of successes.
beta = 1.0;				% Initial scale parameter.
betamin = 1.0e-15; 			% Lower bound on scale.
betamax = 1.0e100;			% Upper bound on scale.
j = 1;					% j counts number of iterations.


% Main optimization loop.
while (j <= options.max_train_iter)

  % Calculate first and second directional derivatives.
  if (success == 1)
    mu = d*gradnew';
    if (mu >= 0)
      d = - gradnew;
      mu = d*gradnew';
    end
    kappa = d*d';
    if kappa < eps
      options(8) = fnow;
      return
    end
    sigma = sigma0/sqrt(kappa);
    xplus = x + sigma*d;
    gplus = strfGrad(strfUnpak(strf,xplus),stim,resp);
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
  fnew = strfErr(strfUnpak(strf,xnew),stim,resp);
  
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

  if options.disp < 0
    if j==1 || mod(j,options.disp)==0
    
        delete([h1 h2]);
        subplot(1,2,1); hold on;
        h1 = bar(x);
        subplot(1,2,2); hold on;
        h2 = plot(options.diagnostics.errs,'.-');
        drawnow;

    end
  elseif options.disp > 0
    if j==1 || mod(j,options.disp)==0
    
      fprintf('iter: %03d | error: %.3f \n',j,fnow);

    end
  end

  if (success == 1)
    % Test for termination

    if (max(abs(alpha*d)) < options.min_step_size & max(abs(fnew-fold)) < options.err_change)
      strf=strfUnpak(strf,x);
      return;

    else
      % Update variables for new position
      fold = fnew;
      gradold = gradnew;
      gradnew = strfGrad(strfUnpak(strf,x),stim,resp);;
     
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
end

% If we get here, then we haven't terminated in the given number of 
% iterations.

strf=strfUnpak(strf,x);
