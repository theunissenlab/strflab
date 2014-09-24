function [strf,err,varargout]=linErr(strf,datIdx)
%function [strf, err, varargout] = linErr(strf, datIdx)
% Evaluate error function for generalized linear model.
%
% INPUTS:
%	  [strf] = a strf stucture
%   [datIdx] = a set of indices into the global [stim] matrix.  Picks out samples to be propagated through model.
%
%
% OUTPUTS:
%      [strf] = strf stucture with updated "internal" field. 
%	    [err] = value of the error function at W for a given STIM and RESP
% [varargout] = (optional) 1xN vector of response before output nonlinearity
%
%	SEE ALSO
%	lin, linPak, linUnPak, linFwd, linGrad
%(Some code modified from NETLAB)
global globDat;

[strf,resp_strf,a]=linFwd(strf,datIdx);
if (size(resp_strf,1) == size(globDat.resp(datIdx),2)) & (size(resp_strf,2) == size(globDat.resp(datIdx),1))
    %This means the resps need to be transposed, that's all.
    resp_strf = resp_strf';
end
switch strf.outputNL

  case 'linear'  	% Linear outputs
      try
    err = 0.5*nansum((resp_strf - globDat.resp(datIdx)).^2);
      catch
          disp(['Error in linErr: possible mismatch in predicted and actual data lengths?']);
          disp(['Size of the STRF response is ' num2str(size(resp_strf)) ' whereas size of real response is ' num2str(size(globDat.resp(datIdx))) '.']);
          error(lasterr);
      end

  case 'logistic'  	% Logistic outputs
    err = - nansum(globDat.resp(datIdx).*log(resp_strf) + (1 - globDat.resp(datIdx)).*log(1 - resp_strf));
   
  case 'softmax'   	% Softmax outputs
    err = - nansum(globDat.resp(datIdx).*log(resp_strf));

  case 'exponential'    % Exponential outputs (for Poisson)
    err =nansum(resp_strf-globDat.resp(datIdx).*log(resp_strf));

  otherwise
    error(['Unknown output nonlinearity ', strf.outputNL]);
end
