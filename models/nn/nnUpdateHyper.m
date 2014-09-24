function strf=nnUpdateHyper(strf,stim,resp)
%NNUPDATEHYPER Re-estimate hyperparameters using evidence approximation.
%
%       Description
%       strf = nnUpdateHyper(strf) re-estimates the hyperparameters ALPHA
%       by applying a Bayesian re-estimation formulae.  The hyperparameter ALPHA can 
%	be a simple scalar if STRF.regtype = 'ridge'
%
%	or a STRF.nin by 1 by # of delays if STRF.regtype='ARD'
%       
%
%	STRF = 	NNUPDATEHYPER(STRF,STIM,RESP) will also update STRF.beta, if the
%		output function is 'linear'.  BETA determines how much to focus on
%		fitting the data compared to the regularization (which STRF.alpha
%		controls)
%
%Some code taken from Netlab

if isfield(strf,'regtype')

   if strcmp(strf.regtype,'ridge')

        w=nnPak(strf);
        strf.alpha=1/mean(w.^2);

   elseif strcmp(strf.regtype,'ARD')

        w1=strf.w1;
        strf.alpha=1./mean(w1.^2,2);
	   
        
   end

   if exist('stim','var')&strcmp(strf.outfn,'linear')
	   reside=resp-nnFwd(strf,stim);
	   strf.beta=1/nanvar(reside);
   end
end
