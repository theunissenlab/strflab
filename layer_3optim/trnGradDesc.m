function [strf,options]=trnGradDescAdaptive(strf,datIdx,options,varargin)

%function [strf,options]=trnGradDescAdaptive(strf,datIdx,options,varargin)
%
% Gradient descent / coordinate descent optimization of STRF.
%
% INPUT:
%                [strf] = model structure obtained via upper level *Init functions
%              [datIdx] = a vector containing indices of the samples to be used in the 
%                         fitting.
%             [options] = option structure containing fields:
%             .funcName = char array, name of this function.
%             .coorDesc = 0 means gradient descent (default)
%                         1 means boosting
%            .earlyStop = 0 means do not early stop (default)
%                         1 means early stop
%              .display = Number of iterations between displaying of errors as we train.
%                         Positive integer means display in stdout.
%                         Negative integer means display in a figure window.
%                         0 means do not display errors (default)
%              .maxIter = Maximum number of training iterations (default: 10000)
%                         if .nIter is supplied, this option has no effect.
%                .nIter = Number of training iterations to use. (default: -1)
%                         Positive integer means .maxIter is ignored, and we fit for this
%                         many iterations, and we return the kernel after the last iteration 
%                         (regardless of the error on the estimation set).  Note that .nIter
%                         can be a nonnegative integer even if options.earlyStop is 1.  In 
%                         this case, we still fit for .nIter iterations, and do not early stop.
%                         -1 means default behavior, ie. stop at maxIter or early stop criterion.
%             .errLastN = Early stopping criterion (default: 30)
%                         +N means stop if we see a series of N iterations which do not
%                         improve performance on the estimation set or do not
%                         improve performance on the stopping set (if there is one).
%                         -N means stop if we see a series of N iterations over which
%                         the slope of a line fit to the estimation error is greater
%                         than options.errSlope OR over which the slope of a line fit
%                         to the stopping error (if there is a stopping set) is greater
%                         than options.errSlope.
%             .errSlope = The slope used to check for convergence (default: -0.001)
%                         (see options.errLastN).
%                         This matters only if options.errLastN is negative.
%             .stepSize = For gradient descent, this is the scale factor for unit-length 
%                         normalized gradient. For boosting, this is the step size.
%                         (default: 0.001)
%             .adaptive = Flag to determine if adaptive step sizes should be used.
%                         (default: 1)
%             .momentum = The non-negative scale factor for previous gradient (default: 0)
%                         Applies only to the gradient descent case.
%                         0 is equivalent to no momentum.
%                         Note that [grad] = N([current grad] + [momentum]*[old grad]), 
%                         where N is unit-length normalization.
%        .nDispTruncate = The number of errors to NOT display at the start.  (Many
%                         errors are swamped until the step size becomes reasonable.
%                         Default = 25.
%                 .is1D = Despaly option to shot the space by time STRF instead of a
%                         bar of the coefficients.  Default 0.
%             .gradNorm = Whether to normalize the gradient at each step (default: 1)
%
%            [varargin] = a vector containing indices of the samples to be used in the 
%	                      early stopping set.
%
%
% OUTPUT:
%                [strf] = structure containing model fit by trnGradDesc.
%
%             [options] = option structure, with additional fields:
%          .diagnostics = contains information on the fitting process.
%
%     .diagnostics.errs = if options.earlyStop is 0,
%                         this has dimensions 1 x N where this has the error on the
%                         training set, over the course of the fitting iterations.
%
%                         if options.earlyStop is 1, this has dimensions 2 x N
%                         where the first row has the error on the training set and
%                         the second row has the error on the stopping set, over
%                         the course of the fitting iterations.
%
% .diagnostics.bestiter = this is the iteration number (1-indexed) corresponding to the
%                         returned kernel.
%
% EXAMPLE:
% options = trnGradDesc
%   returns default options
%
% SEE ALSO: trnSCG, trnDirectFit
%
%(Some code modified from NETLAB)
%
% nDispTruncate Bug Fix, adaptive off by default, add titles to graphs - by Michael Oliver (June 2009)
% ====================


% Set default option values
% --------------------
global globDat
optDeflt.funcName='trnGradDesc';
optRange.funcName={'trnGradDesc'};

optDeflt.display=0;
optRange.display=[-1e6 1e6];
optDeflt.coorDesc=0;
optRange.coorDesc=[0 1];
optDeflt.earlyStop=double([length(varargin)>0]);
optRange.earlyStop=[0 1];
optDeflt.stepSize=0.00001;
optRange.stepSize=double([eps('single') Inf]);
optDeflt.maxIter=10000;
optRange.maxIter=[0 Inf];
optDeflt.nIter=-1;
optRange.nIter=[-1 Inf];
optDeflt.errLastN=30;
optRange.errLastN=[-Inf Inf];
optDeflt.errSlope=-0.001;
optRange.errSlope=[-Inf Inf];
optDeflt.momentum=0;
optRange.momentum=[0 Inf];
optDeflt.adaptive = 0;
optRange.adaptive = [0 1];
optDeflt.nDispTruncate = 25;
optRange.nDispTruncate = [0 Inf];
optDeflt.is1D = 0;
optRange.is1D = [0 1];
optDeflt.gradNorm = 1;
optRange.gradNorm = [0 1];
if nargin<3
  options=optDeflt;
else
  options=defaultOpt(options,optDeflt,optRange);
end
if nargin<1
  strf=optDeflt;
  return;
end

if options.earlyStop & isempty(varargin)
    warning(sprintf('You must specify an Early Stopping Index if you want to use Early Stopping'));
    return
end

% calc
if options.nIter == -1
  iterstouse=options.maxIter;
else
  iterstouse=options.nIter;
end


% initialize plot
if options.display < 0
    clf;
    if options.earlyStop==0
        set(gcf,'Units','points','Position',[20 50 800 300]);
        subplot(1,2,1);
        if options.is1D == 0
            xlabel('parameter');
            ylabel('magnitude');
        end
        subplot(1,2,2);
        xlabel('iteration');
        ylabel('estimation error');
    else
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
end


% init
options.diagnostics.errs = [];  % record of estimation set error
esterr_min = Inf;               % minimum estimation error found so far
stoperr_min = Inf;              % minimum stopping error found so far
grad_prev = 0;                  % previous gradient
h1 = []; h2 = []; h3 = [];      % plot handles
[strf,x] = strfPak(strf);              % initial seed

estbadcnt = 0;
stopbadcnt = 0;

%Initialize adaptive gradient parameters

flagValue = -99999;
thisErr = flagValue;
lastErr = flagValue;
n_decreasing = 1;
% do the loop
for iter=1:iterstouse
  
  % if not the first iteration, adjust according to the gradient
  if iter~=1

    % calculate gradient
    [strf,grad] = strfGrad(strf,datIdx);

    % adjust the kernel
    if options.coorDesc==0
      if options.gradNorm
        grad = grad/norm(grad);
      end
      if options.momentum
        grad = grad + options.momentum * grad_prev;
        grad = grad/norm(grad);
      end
  
      x = x - (options.stepSize * grad);
    else
      [d,idx] = max(abs(grad));
      x(idx) = x(idx) - options.stepSize * sign(grad(idx));
    end

    % record gradient for next iteration
    grad_prev = grad;
    
  end
  
  % calculate and record error
  strf=strfUnpak(strf,x);
  [strf,esterr] = strfErr(strf,datIdx);
  if options.adaptive
      thisErr = esterr;
      if any([thisErr lastErr] == flagValue)
          % Can't tell if improving yet
      else
          if lastErr < thisErr
              %Step size too large.  Decrease it.
              %Optimization is fastest if the step size is roughly half of the
              %smallest step size which results in an error increase, so it's
              %fastest to hang around half of this first-smallest-bad step size
              %for as long as possible.
              options.stepSize = options.stepSize *(.5^n_decreasing);
              n_decreasing = n_decreasing + 1;
          else
              %Step size could be bigger
              options.stepSize = options.stepSize*(1.1);  %This term will grow slowly after it's been last cut in half, but then quickly find a too-big step.
              n_decreasing = 1;
          end
          %disp(['New step size is ' num2str(options.stepSize) '.']);
      end
      lastErr = thisErr;
  end
  options.diagnostics.errs(1,iter) = esterr;
  if options.earlyStop
    [strf,stoperr] = strfErr(strf,varargin{:});
    options.diagnostics.errs(2,iter) = stoperr;
  end

  % report
  if options.display > 0
    if iter==1 || mod(iter,options.display)==0
      if options.earlyStop==0
        fprintf('iter: %03d | est: %.3f (%s)\n',iter,esterr,choose(esterr < esterr_min,'D','U'));
      else
        fprintf('iter: %03d | est: %.3f (%s) | stop: %.3f (%s)\n',iter,esterr,choose(esterr < esterr_min,'D','U'), ...
          stoperr,choose(stoperr < stoperr_min,'D','U'));
      end
    end
  end
  if options.display < 0
    if iter==1 || mod(iter,options.display)==0
      if options.earlyStop==0
        delete([h1 h2]);
        subplot(1,2,1); hold on; %axis tight;
        if options.is1D
            simagesc(reshape(x(1:(end-1)),strf.nIn,length(strf.delays)));
            axis tight;
        else
            
            h1 = bar(x(1:end-1)); title('Parameter magnitudes');
        end
        subplot(1,2,2); hold on;
        h2 = plot(options.diagnostics.errs,'-'); title('Training Set Error');
      else
        delete([h1; h2; h3]);
        subplot(1,3,1); hold on; %axis tight;
        if options.is1D
            simagesc(reshape(x(1:(end-1)),strf.nIn,length(strf.delays)));
            axis tight;
        else
            h1 = bar(x(1:end-1)); title('Parameter magnitudes');
        end
        subplot(1,3,2); hold on;
        h2 = plot((1+options.nDispTruncate):length(options.diagnostics.errs),options.diagnostics.errs(1,(1+options.nDispTruncate):end),'-'); title('Training Set Error');
        subplot(1,3,3); hold on;
        h3 = plot((1+options.nDispTruncate):length(options.diagnostics.errs),options.diagnostics.errs(2,(1+options.nDispTruncate):end),'-'); title('Early Stopping Set Error');
      end
      drawnow;
    end
  end

  % do we consider this iteration to be the best yet?
  % (if fixed iter is supplied) OR (if no early stopping and estimation error is minimum so far) OR (if early stopping and stopping error is minimum so far)
  if (options.nIter >= 0) || (options.earlyStop==0 && esterr < esterr_min) || (options.earlyStop==1 && stoperr < stoperr_min)
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
  if options.earlyStop
    if stoperr < stoperr_min
      stopbadcnt = 0;
      stoperr_min = stoperr;
    else
      stopbadcnt = stopbadcnt + 1;
    end
  end
  
  % do we stop? check the slopes.
  % (if not doing the special fixed iter case) AND (if want the slope check) && (we have enough iterations)
  if (options.nIter == -1) && (options.errLastN < 0) && (iter - (-options.errLastN) + 1 >= 1)

    % check estimation error slope (both cases)
    params = polyfit(1:-options.errLastN,options.diagnostics.errs(1,iter-(-options.errLastN)+1 : end),1);
    if params(1) > options.errSlope
      break;
    end
    
    % check stopping error slope (only for early stopping)
    if options.earlyStop
      params = polyfit(1:-options.errLastN,options.diagnostics.errs(2,iter-(-options.errLastN)+1 : end),1);
      if params(1) > options.errSlope
        break;
      end
    end
    
  end
  
  % do we stop? check the number of iterations since improvement
  % (if not doing the special fixed iter case) AND (either the estimation error or stopping error has bottomed out)
  if (options.nIter == -1) && (options.errLastN > 0 && (estbadcnt == options.errLastN || stopbadcnt == options.errLastN))
    break;
  end

end

% output
strf = strfUnpak(strf,x_best);

