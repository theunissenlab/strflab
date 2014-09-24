function [strf, options] = trnLARS(strf, datIdx, options, varargin)

% function [strf, options] = trnLARS(strf, datIdx, options, varargin)
%
% Least angle regression optimization for STRF. Stimulus matrix should be z-scored
%
% INPUT:
%                [strf] = model structure obtained via upper level *Init functions
%              [datIdx] = a vector containing indices of the samples to be used in the 
%                         fitting.
%             [options] = option structure containing fields:
%             .funcName = char array, name of this function.
%               .method = 'lars' for least angle regression (default).
%                         'en' for elastic net regression. 
%                         'lasso' for the LASSO modification.
%                 .type = 'linear' for regular lars, in case STRF is in 'glm' class.
%                         'product' for STRFs in 'oil' class. Includes products 
%						  between the regressors in the stimulus matrix. The algorithm 
%					      iteratively includes multiplicative terms as variables join
%                         the active set (default) 
%                         'rectproduct' for STRFs in 'oil' class. Creates 4 parameters 
%						  for each pair of non-linear combinations: separates each vector 
%                         into positive and negative parts and makes all possible
%                   	  combinations between them [(+,+);(-,-);(+,-);(-,+)]
%              .lambda2 = User-defined ridge term coefficient, for elastic net regression
%						  0 means elastic net will not be performed (default)
%                .maxL1 = Upper bound on L1 norm of the regression coefficients for lasso
%						  0 means lasso will not be performed (default)
%            .earlyStop = 0 means do not early stop 
%                         1 means early stop (default)
%                         -1 means early stop by calculating estimation risk    
%              .display = Number of iterations between displaying of errors as we train.
%                         Positive integer means display in stdout.
%                         Negative integer means display in a figure window.
%                         0 means do not display errors (default)
%              .maxIter = Maximum number of training iterations (default: 10000)
%                         if .nIter is supplied, this option has no effect.
%              .maxChan = Maximum number of variables allowed to enter the active set
%						  0 means no bounds (default)
%                         > 0 indicates how many variables
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
%        .nDispTruncate = The number of errors to NOT display at the start.  (Many
%                         errors are swamped until the step size becomes reasonable.
%                         Default = 25.
%                 .is1D = Despaly option to shot the space by time STRF instead of a
%                         bar of the coefficients.  Default 0.
%            [varargin] = a vector containing indices of the samples to be used in the 
%	                      early stopping set.
%
%
% OUTPUT:
%                [strf] = structure containing model fit by trnGradDesc.
%            .varlookup = variable lookup table where first column
%                         corresponds to stimulus column and second column
%                         corresponds to delay
%           .prod_pairs = lookup table for non-linear interaction terms
%                         if type is 'product', this is a p by 2 matrix where
%                         p is the number of non-linear terms added to stim.
%                         Each row contains the indices of the terms from the
%                         original stim that where multiplied, in the order in
%                         which they enter the active set.
%                         if type is 'rectproduct', this is a p by 4 matrix, 
%                         where p is as above. In addition to the indices of 
%                         the terms from the original X, the two other elements 
%                         in each row are 1 or -1 if the equivalent parameter
%                         vector has just positive or negative values,
%                         respectively.
%          .interaction = for STRFs of the 'oil' type, string specifying type of 
%                         non-linear interaction
%
%                   [options] = option structure, with additional fields:
%                .diagnostics = contains information on the fitting process.
%
%           .diagnostics.errs = if options.earlyStop is 0,
%                               this has dimensions 1 x N where this has the error on the
%                               training set, over the course of the fitting iterations.
%
%                               if options.earlyStop is 1, this has dimensions 2 x N
%                               where the first row has the error on the training set and
%                               the second row has the error on the stopping set, over
%                               the course of the fitting iterations.
%
%       .diagnostics.bestiter = this is the iteration number (1-indexed) corresponding to the
%                               returned kernel.
%  .diagnostics.nonZeroWtsIdx = 2 column matrix where the first contains the
%                               indices of non-zero weights in beta and the second
%                               column contains the respective weight
%                               values
%         .diagnostics.linIdx = non-zero weights of linear terms, as above
%         .diagnostics.linNum = number of non-zero linear term weights
%        .diagnostics.nlinIdx = non-zero weights of non-linear terms, as above
%        .diagnostics.nlinNum = number of non-zero non-linear term weights
%
% SEE ALSO: computeLARSstep, RchangeChol, earlyStopTrim, changeVarIdx, oilAddStim
%
% Modified from the implementation of Karl Skoglund, IMM, DTU, kas@imm.dtu.dk
% (@ http://www2.imm.dtu.dk/pubdb/views/publication_details.php?id=3897)
% References: 1. 'Least Angle Regression' by B. Efron et al, 2004. Ann Stats
% 2. 'Regularization and Variable Selection via the Elastic Net' by H. Zou & T. Hastie, 2005.

%% Input checking
global globDat
optDeflt.funcName='trnLARS';
optRange.funcName={'trnLARS'};

optDeflt.method= 'lars';
optRange.method={'lars', 'en', 'lasso'};
optDeflt.type= 'linear';
optRange.type={'product','inv_product','pos_product','rectproduct', 'linear', ...
 	'division', 'divnorm', 'sqrt_sqsum', 'sqsum', 'softmax' ,'all'};
optDeflt.lambda2 = 0;
optRange.lambda2 = [0 Inf];
optDeflt.maxL1=0;
optRange.maxL1=[0 Inf];
optDeflt.earlyStop=1;
optRange.earlyStop=[-1 1];
optDeflt.display=1;
optRange.display=[-1e6 1e6];
optDeflt.maxIter=500;
optRange.maxIter=[0 Inf];
optDeflt.maxChan=0;
optRange.maxChan=[0 Inf];
optDeflt.nIter=-1;
optRange.nIter=[-1 Inf];
optDeflt.errLastN=30;
optRange.errLastN=[-Inf Inf];
optDeflt.errSlope=-0.001;
optRange.errSlope=[-Inf Inf];
optDeflt.nDispTruncate = 25;
optRange.nDispTruncate = [0 Inf];
optDeflt.is1D = 0;
optRange.is1D = [0 1];

if nargin<3
  options=optDeflt;
else
  options=defaultOpt(options,optDeflt,optRange);
end
if nargin<1
  strf=optDeflt;
  return;
end

lambda2 = options.lambda2;

%% Variable setup
n = length(datIdx); p = strf.nIn;

if strcmpi(strf.type,'glm')
	options.type = 'linear';
end

strf.prod_pairs = []; % variable containing the index of pairs whose product is in stim

strf.varlookup = [];
for iii = 1:length(strf.delays)
    strf.varlookup=[strf.varlookup; (1:strf.nIn)' ones(size((1:strf.nIn)'))*strf.delays(iii)];
end

if ~strcmpi(options.type,'linear') && strcmpi(strf.type,'oil')
    strf.interaction = options.type;
else strf.interaction = [];
end

if strcmpi(options.type,'divnorm') || strcmpi(options.type,'all') 
    stimsum = sum(globDat.stim,2);
end

options.diagnostics.errs = [];  % record of estimation set error
options.diagnostics.risk = [];
options.diagnostics.betas = [];	% cell array containing the beta vector for each iteration (k+1)
esterr_min = Inf;               % minimum estimation error found so far
stoperr_min = Inf;              % minimum stopping error found so far
risk_min = Inf;
h1 = []; h2 = []; h3 = [];      % plot handles

estbadcnt = 0;
stopbadcnt = 0;
riskbadcnt = 0;

ndelays = length(strf.delays);
nvars = (strf.nWts-1);

if options.earlyStop && options.maxL1 > 0,
  options.maxL1 = options.maxL1/sqrt(1 + lambda2);
end

if ~options.earlyStop
  beta = zeros(options.maxIter,nvars);
elseif options.earlyStop && options.maxChan > 0
  beta = zeros(2*round(options.maxChan), nvars);
else
  beta = zeros(200, nvars);
end

mu = zeros(n, 1); % current "position" as LARS-EN travels towards lsq solution
I = 1:p*ndelays; % inactive set
A = []; % active set
R = []; % Cholesky factorization R'R = X'X where R is upper triangular

lassocond = 0; % Set to 1 if LASSO condition is met
stopcond = 0; % Set to 1 if early stopping condition is met
iter = 0; % Algorithm step count
vars = 0; % Current number of variables

if strcmp(options.method, 'en') && options.lambda2 == 0
    error('You should set lambda2 > 0 for elastic net')
end

if ~strcmp(options.method, 'en') && options.lambda2 ~= 0
    error('You should set lambda2 = 0 for lasso and lars')
end

if strcmp(options.method, 'en')
    d1 = sqrt(lambda2); % Convenience variables d1 and d2
    d2 = 1/sqrt(1 + lambda2); 
else
    d1 = 0;
    d2 = 1;
end

% Calculate baseline to be subtracted so we work with centered responses
if strf.b1 == 0
    strf.b1 = mean(globDat.resp(datIdx));
	globDat.resp(datIdx) = globDat.resp(datIdx) - strf.b1;
end

%% Added
strf.orig_nIn = strf.nIn;

%% Main loop
while vars < nvars && ~stopcond && iter < options.maxIter

     if ~strcmpi(options.type,'linear')
		 if strcmpi(options.type,'divnorm') && vars >= 1
			[strf, beta, A, I] = oilAddStim(strf,beta,options,A,I,iter,stimsum);
		 elseif strcmpi(options.type,'all') && vars >= 2
			[strf, beta, A, I] = oilAddStim(strf,beta,options,A,I,iter,stimsum);
         elseif ~strcmpi(options.type,'divnorm') && vars >= 2 % if there are at least 2 variables in the active set, include their products in stim and update   
             [strf, beta, A, I] = oilAddStim(strf,beta,options,A,I,iter);
         end
     end
    
  iter = iter + 1;
  [strf, grad] = strfGrad(strf,datIdx);
  grad = -grad(1:end-1)*d2;
  [C j] = max(abs(grad(I)));
  j = I(j);
  if ~lassocond % if a variable has been dropped, do one iteration with this configuration (don't add new one right away)
    R = RchangeChol(R,'insert',j,strf,A,lambda2,datIdx);
    A = [A j];
    I(I == j) = [];
    vars = vars + 1;
  else
    j=0;
  end

  [gamma,u1,w] = computeLARSstep(grad,R,A,I,C,d1,d2,datIdx,strf); 

  if ~strcmp(options.method, 'lars')
      % LASSO modification
      lassocond = 0;
      temp = -beta(iter,A)./w';
      [gamma_tilde] = min([temp(temp > 0) gamma]);
      j2 = find(temp == gamma_tilde);
      if gamma_tilde < gamma,
        gamma = gamma_tilde;
        lassocond = 1;
      end
  end

  mu = mu + gamma*u1;
  if size(beta,1) < iter+1
    beta = [beta; zeros(size(beta))];
  end
  % figure out how to go from A to beta
  beta(iter+1,A) = beta(iter,A) + gamma*w'; 
  strf = strfUnpak(strf,[beta(iter+1,:) strf.b1]);

  % calculate and record error
  % strf=strfUnpak(strf,x);
  [strf,esterr] = strfErr(strf,datIdx);
    
  options.diagnostics.errs(1,iter) = esterr;
  options.diagnostics.betas{iter} = beta(iter+1,:);
    %%added == 1 so error not checked when its -1 
  if options.earlyStop == 1
      [strf,stoperr] = strfErr(strf,varargin{:});
      options.diagnostics.errs(2,iter) = stoperr;
  end
  if options.earlyStop == -1
      Cp = LARSrisk(esterr,vars,datIdx);
      options.diagnostics.risk(1,iter) = Cp;
  end
  
  % If LASSO condition satisfied, drop variable from active set
  if lassocond == 1
    R = RchangeChol(R,'delete',j2); %R = choldelete(R,j);
    I = [I A(j2)];
    jd = A(j2);
    A(j2) = [];
    vars = vars - 1;
  else
    jd = 0;
  end
  
  if options.display > 0
      if iter==1 || mod(iter,options.display) == 0
          if options.earlyStop == 0
              fprintf('iter: %04d | added: %05d | dropped: %05d | active set: %03d | est: %.3f (%s) | stepsize: %.3f \n',...
                  iter, j, jd, vars, esterr,choose(esterr < esterr_min,'D','U'), gamma);
          elseif options.earlyStop == 1
              fprintf('iter: %04d | added: %05d | dropped: %05d | active set: %03d | est: %.3f (%s) | stop: %.3f (%s) | stepsize: %.3f \n',...
                  iter, j, jd, vars, esterr, choose(esterr < esterr_min,'D','U'),stoperr,choose(stoperr < stoperr_min,'D','U'), gamma);
          elseif options.earlyStop == -1
              fprintf('iter: %04d | added: %05d | dropped: %05d | active set: %03d | est: %.3f (%s) | risk: %.3f (%s) | stepsize: %.3f \n',...
                  iter, j, jd, vars, esterr, choose(esterr < esterr_min,'D','U'),Cp,choose(Cp < risk_min,'D','U'), gamma);
          end
      end
  end
  
% report
  if options.display < 0
    if iter==1 || mod(iter,options.display)==0
      if options.earlyStop==0
        delete([h1 h2]);
        subplot(1,2,1); hold on; %axis tight;
        if options.is1D
            simagesc(reshape(beta(iter+1,:),strf.nIn,length(strf.delays)));
            axis tight;
        else         
            h1 = bar(beta(iter+1,:)); title('Parameter magnitudes');
        end
        subplot(1,2,2); hold on;
        h2 = plot(options.diagnostics.errs,'-'); title('Training Set Error');
      else
        delete([h1; h2; h3]);
        subplot(1,3,1); hold on; %axis tight;
        if options.is1D
            simagesc(reshape(beta(iter+1,:),strf.nIn,length(strf.delays)));
            axis tight;
        else
            h1 = bar(beta(iter+1,:)); title('Parameter magnitudes');
        end
        subplot(1,3,2); hold on;
        h2 = plot((1+options.nDispTruncate):length(options.diagnostics.errs),...
            options.diagnostics.errs(1,(1+options.nDispTruncate):end),'-'); title('Training Set Error');
        subplot(1,3,3); hold on;
        h3 = plot((1+options.nDispTruncate):length(options.diagnostics.errs),...
            options.diagnostics.errs(2,(1+options.nDispTruncate):end),'-'); title('Early Stopping Set Error');
      end
      drawnow;
    end
  end

  % do we consider this iteration to be the best yet?
  % (if fixed iter is supplied) OR (if no early stopping and estimation error is minimum so far) OR (if early stopping and stopping error is minimum so far)
  if (options.nIter >= 0) || (options.earlyStop==0 && esterr < esterr_min) || ...
          (options.earlyStop==1 && stoperr < stoperr_min) || (options.earlyStop==-1 && Cp < risk_min)
    options.diagnostics.bestiter = iter;
    options.diagnostics.bestiter_varlookup = strf.varlookup;
    options.diagnostics.bestiter_prod_pairs = strf.prod_pairs;
    options.diagnostics.bestbeta = beta;
  end

  % check estimation error against minimum (both cases)
  if esterr < esterr_min
    estbadcnt = 0;
    esterr_min = esterr;
  else
    estbadcnt = estbadcnt + 1;
  end

  % check stopping error against minimum (only for early stopping)
  if options.earlyStop == 1
    if stoperr < stoperr_min
      stopbadcnt = 0;
      stoperr_min = stoperr;
    else
      stopbadcnt = stopbadcnt + 1;
    end
  end
  if options.earlyStop == -1
     if Cp < risk_min
      riskbadcnt = 0;
      risk_min = Cp;
    else
      riskbadcnt = riskbadcnt + 1;
    end     
  end
  
  % do we stop? check the slopes.
  % (if not doing the special fixed iter case) AND (if want the slope check) && (we have enough iterations)
  if (options.nIter == -1) && (options.errLastN < 0) && (iter - (-options.errLastN) + 1 >= 1)

    % check estimation error slope (both cases)
    params = polyfit(1:-options.errLastN,options.diagnostics.errs(1,iter-(-options.errLastN)+1 : end),1);
    if params(1) > options.errSlope
      stopcond = 1;
    end
    
    % check stopping error slope (only for early stopping)
    if options.earlyStop || options.earlyStop == -1
        if options.earlyStop
            params = polyfit(1:-options.errLastN,options.diagnostics.errs(2,iter-(-options.errLastN)+1 : end),1);
        elseif options.earlyStop == -1
            params = polyfit(1:-options.errLastN,options.diagnostics.risk(1,iter-(-options.errLastN)+1 : end),1);
        end
        if params(1) > options.errSlope
            stopcond = 1;
        end
    end
    
  end
  
  % do we stop? check the number of iterations since improvement
  % (if not doing the special fixed iter case) AND (either the estimation error or stopping error has bottomed out)
  if (options.nIter == -1) && (options.errLastN > 0 && (estbadcnt == options.errLastN || stopbadcnt == options.errLastN))
    stopcond = 1;
  end
  if options.earlyStop == -1 && (options.nIter == -1) && (options.errLastN > 0 && (riskbadcnt == options.errLastN))
      stopcond = 1;
  end
  
  if stopcond == 1
      if ~strcmpi(options.type,'linear') % delete non-linear terms included over last .errLastN iterations
            [strf,beta] = earlyStopTrim(strf,options);
      end
  end

  % Early stopping at specified bound on L1 norm of beta
  if options.earlyStop & options.maxL1 > 0
    t2 = sum(abs(beta(iter+1,:)));
    if t2 >= options.maxL1
      t1 = sum(abs(beta(iter,:)));
      s = (options.maxL1 - t1)/(t2 - t1); % interpolation factor 0 < s < 1
      beta(iter+1,:) = beta(iter,:) + s*(beta(iter+1,:) - beta(iter,:));
      stopcond = 1;
    end
  end

  % Early stopping at specified number of variables
  if options.earlyStop & options.maxChan > 0
    stopcond = vars >= options.maxChan;
  end

end

% trim beta 
if size(beta,1) > iter+1
  beta(iter+2:end, :) = [];
end

if strcmp(options.method, 'en')
    beta = beta/d2;    % divide by d2 to avoid double shrinkage
end

% output
strf = strfUnpak(strf,[beta(options.diagnostics.bestiter+1,:) strf.b1]);

options.diagnostics.nonZeroWtsIdx = (find(options.diagnostics.bestbeta(options.diagnostics.bestiter+1,:) ~= 0))';
options.diagnostics.nonZeroWtsIdx = [options.diagnostics.nonZeroWtsIdx ...
    (options.diagnostics.bestbeta(options.diagnostics.bestiter+1,options.diagnostics.nonZeroWtsIdx))'];
options.diagnostics.linIdx = find(strf.varlookup(options.diagnostics.nonZeroWtsIdx(:,1),1) <= strf.orig_nIn);
options.diagnostics.linIdx = [options.diagnostics.linIdx  options.diagnostics.nonZeroWtsIdx(options.diagnostics.linIdx,2)]; 
options.diagnostics.linNum = size(options.diagnostics.linIdx,1);
options.diagnostics.nlinIdx = find(strf.varlookup(options.diagnostics.nonZeroWtsIdx(:,1),1) > strf.orig_nIn);
options.diagnostics.nlinIdx = [options.diagnostics.nlinIdx  options.diagnostics.nonZeroWtsIdx(options.diagnostics.nlinIdx,2)]; 
options.diagnostics.nlinNum = size(options.diagnostics.nlinIdx,1);

if iter == options.maxIter
  disp('LARS-EN warning: Forced exit. Maximum number of iteration reached.');
end