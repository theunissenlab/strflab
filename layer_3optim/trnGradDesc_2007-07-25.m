function [strf, options] = gradDescTrain(strf, stim, resp, options, varargin)

%GRADDESCTRAIN Gradient descent optimization.
%
% [STRF, OPTIONS] = GRADDESCTRAIN(STRF, STIM, RESP)
% [STRF, OPTIONS] = GRADDESCTRAIN(STRF, STIM, RESP, OPTIONS)
% [STRF, OPTIONS] = GRADDESCTRAIN(STRF, STIM, RESP, OPTIONS, STIMVAL, RESPVAL)
%        OPTIONS  = GRADDESCTRAIN
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
% OPTIONS.descenttype (nonnegative integer)
%
%   0 means gradient descent
%   1 means boosting
%   (default, 0).
%
% OPTIONS.earlystop (nonnegative integer)
%   0 means do not early stop
%   1 means early stop
%   (default, 0).
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
% OPTIONS.fixed_iter (-1 or nonnegative integer)
%   Number of training iterations to use.  If nonnegative integer,
%   .max_train_iter is ignored, and we fit for this many iterations,
%   and we return the kernel after the last iteration (regardless of
%   the error on the estimation set).  Note that .fixed_iter can be
%   a nonnegative integer even if options.earlystop is 1.  In this case,
%   we still fit for .fixed_iter iterations, and do not early stop.
%   If -1, this means default behavior.  (default, -1).
%
% OPTIONS.err_change (positive integer or negative integer)
%   If N, stop if we see a series of N iterations which do not
%     improve performance on the estimation set or do not
%     improve performance on the stopping set (if there is one).
%   If -N, stop if we see a series of N iterations over which
%     the slope of a line fit to the estimation error is greater
%     than options.err_slope or over which the slope of a line fit
%     to the stopping error (if there is a stopping set) is greater
%     than options.err_slope.
%   (default, 30).
%
% OPTIONS.err_slope (number)
%   The slope used to check for convergence (see options.err_change).
%   This matters only if options.err_change is negative.
%   (default, -0.001).
%
% OPTIONS.stepsize (positive number)
%
%   For gradient descent, this is the scale factor for unit-length normalized gradient.
%   For boosting, this is the step size.
%   (default, 0.001).
%
% OPTIONS.momentum (nonnegative number)
%
%   Applies only to the gradient descent case.
%   This is the scale factor for previous gradient.  If 0, this is equivalent to no momentum.
%   Note that [grad] = N([current grad] + [momentum]*[old grad]), where N is unit-length normalization.
%   (default, 0).
%
% STIMVAL (provide if and only if options.earlystop is 1)
%
% RESPVAL (provide if and only if options.earlystop is 1)
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
%   if options.earlystop is 0,
%   this has dimensions 1 x N where this has the error on the
%   training set, over the course of the fitting iterations.
%
%   if options.earlystop is 1, this has dimensions 2 x N
%   where the first row has the error on the training set and
%   the second row has the error on the stopping set, over
%   the course of the fitting iterations.
%
% OPTIONS.diagnostics.bestiter
%
%   this is the iteration number (1-indexed) corresponding to the
%   returned kernel.
% 
%(Some code modified from NETLAB)

% defaults
optDeflt.descenttype = 0;
optRange.descenttype = [0 1];
optDeflt.earlystop = 0;
optRange.earlystop = [0 1];
optDeflt.disp = 0;
optRange.disp = [-Inf Inf];
optDeflt.max_train_iter = 1000;
optRange.max_train_iter = [0 Inf];
optDeflt.fixed_iter = -1;
optRange.fixed_iter = [-1 Inf];
optDeflt.err_change = 30;
optRange.err_change = [-Inf Inf];
optDeflt.err_slope = -0.001;
optRange.err_slope = [-Inf Inf];
optDeflt.stepsize = 0.001;
optRange.stepsize = [eps Inf];
optDeflt.momentum = 0;
optRange.momentum = [0 Inf];
if nargin < 4
  options = optDeflt;
else
  options = defaultOpt(options,optDeflt,optRange);
end
if nargin < 3
  strf = optDeflt;
  return;
end

% calc
if options.fixed_iter == -1
  iterstouse = options.max_train_iter;
else
  iterstouse = options.fixed_iter;
end

% initialize plot
if options.disp < 0
  figure;
  set(gcf,'Units','points','Position',[100 100 800 400]);
  if options.earlystop==0
    subplot(1,2,1);
    xlabel('parameter');
    ylabel('magnitude');
    subplot(1,2,2);
    xlabel('iteration');
    ylabel('estimation error');
  else
    subplot(1,3,1);
    xlabel('parameter');
    ylabel('magnitude');
    subplot(1,3,2);
    xlabel('iteration');
    ylabel('estimation error');
    subplot(1,3,3);
    xlabel('iteration');
    ylabel('stopping error');
  end
end

% init
options.diagnostics.errs = [];  % record of estimation set error
esterr_min = Inf;               % minimum estimation error found so far
stoperr_min = Inf;              % minimum stopping error found so far
grad_prev = 0;                  % previous gradient
h1 = []; h2 = []; h3 = [];      % plot handles
x = strfPak(strf);              % initial seed

estbadcnt = 0;
stopbadcnt = 0;


% do the loop
for iter=1:iterstouse
  
  % if not the first iteration, adjust according to the gradient
  if iter~=1

    % calculate gradient
    grad = strfGrad(strf,stim,resp);

    % adjust the kernel
    if options.descenttype==0
      grad = grad/norm(grad);
      if options.momentum
        grad = grad + options.momentum * grad_prev;
        grad = grad/norm(grad);
      end
  
      x = x - (options.stepsize * grad);
    else
      [d,idx] = max(abs(grad));
      x(idx) = x(idx) - options.stepsize * sign(grad(idx));
    end

    % record gradient for next iteration
    grad_prev = grad;
    
  end
  
  % calculate and record error
  strf=strfUnpak(strf,x);
  esterr = strfErr(strf,stim,resp);
  options.diagnostics.errs(1,iter) = esterr;
  if options.earlystop
    stoperr = strfErr(strf,varargin{:});
    options.diagnostics.errs(2,iter) = stoperr;
  end

  % report
  if options.disp > 0
    if iter==1 || mod(iter,options.disp)==0
      if options.earlystop==0
        fprintf('iter: %03d | est: %.3f (%s)\n',iter,esterr,choose(esterr < esterr_min,'D','U'));
      else
        fprintf('iter: %03d | est: %.3f (%s) | stop: %.3f (%s)\n',iter,esterr,choose(esterr < esterr_min,'D','U'), ...
          stoperr,choose(stoperr < stoperr_min,'D','U'));
      end
    end
  end
  if options.disp < 0
    if iter==1 || mod(iter,options.disp)==0
      if options.earlystop==0
        delete([h1 h2]);
        subplot(1,2,1); hold on;
        h1 = bar(x);
        subplot(1,2,2); hold on;
        h2 = plot(options.diagnostics.errs,'.-');
      else
        delete([h1 h2 h3]);
        subplot(1,3,1); hold on;
        h1 = bar(x);
        subplot(1,3,2); hold on;
        h2 = plot(options.diagnostics.errs(1,:),'.-');
        subplot(1,3,3); hold on;
        h3 = plot(options.diagnostics.errs(2,:),'.-');
      end
      drawnow;
    end
  end

  % do we consider this iteration to be the best yet?
  % (if fixed iter is supplied) OR (if no early stopping and estimation error is minimum so far) OR (if early stopping and stopping error is minimum so far)
  if (options.fixed_iter >= 0) || (options.earlystop==0 && esterr < esterr_min) || (options.earlystop==1 && stoperr < stoperr_min)
    options.diagnostics.bestiter = iter;
    x_best = x;
  end

  % check estimation error against minimum (both cases)
  if esterr < esterr_min
    estbadcnt = 0;
    esterr_min = esterr;
  else
    estbadcnt = estbadcnt + 1;
  end

  % check stopping error against minimum (only for early stopping)
  if options.earlystop
    if stoperr < stoperr_min
      stopbadcnt = 0;
      stoperr_min = stoperr;
    else
      stopbadcnt = stopbadcnt + 1;
    end
  end
  
  % do we stop? check the slopes.
  % (if not doing the special fixed iter case) AND (if want the slope check) && (we have enough iterations)
  if (options.fixed_iter == -1) && (options.err_change < 0) && (iter - (-options.err_change) + 1 >= 1)

    % check estimation error slope (both cases)
    params = polyfit(1:-options.err_change,options.diagnostics.errs(1,iter-(-options.err_change)+1 : end),1);
    if params(1) > options.err_slope
      break;
    end
    
    % check stopping error slope (only for early stopping)
    if options.earlystop
      params = polyfit(1:-options.err_change,options.diagnostics.errs(2,iter-(-options.err_change)+1 : end),1);
      if params(1) > options.err_slope
        break;
      end
    end
    
  end
  
  % do we stop? check the number of iterations since improvement
  % (if not doing the special fixed iter case) AND (either the estimation error or stopping error has bottomed out)
  if (options.fixed_iter == -1) && (options.err_change > 0 && (estbadcnt == options.err_change || stopbadcnt == options.err_change))
    break;
  end

end

% output
strf = strfUnpak(strf,x_best);
