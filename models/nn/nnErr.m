function [strf,err,varargout] = nnErr(strf,datIdx)
%NNERR	Evaluate error function for 2-layer network.
%
%	Description
%	E = NNERR(STRF, STIM, RESP) takes a network data structure STRF 
%	together with a matrix X of input vectors and a matrix T of target 
%	vectors, and evaluates the error function E. The choice of error 
%	function corresponds to the output unit activation function. Each 
%	row of X corresponds to one input vector and each row of T 
%	corresponds to one target vector.
%
%	[E, VARARGOUT] = NNERR(STRF, STIM, RESP) 
%
%Some code taken from Netlab

global globDat;

[strf,resp_strf]=nnFwd(strf,datIdx);

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


[err, eprior] = nnErrBayes(strf, err);
