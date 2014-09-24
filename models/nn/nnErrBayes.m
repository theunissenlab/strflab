function [e, eprior]=nnErrBayes(strf, e)
%NNERR	Evaluate error function for 2-layer network.
%
%	Description
%	[E, EPRIOR]= NNERRBAYES(STRF, E) takes a network data structure STRF
%	and an error term E, and adds the error due to regularization.
%	
%
%	Inputs:
%	
%	STRF, a neural network strf structure
%		nnErrBayes uses the field strf.regtype to determine how
%		which method of regularization to use. 
%
%		strf.regtype='ridge' penalizes the squared sum of all the
%			             strf weights times strf.alpha
%
%		strf.regtype='ARD' penalizes the weights of each input
%				   channel independently
%
%		strf.alpha, the hyperparameters which determines how much to
%			    penalize the weights.
%		
%		strf.beta,  the hyperparameter that dtermines how much to
%			    penalize misfitting the data.
%
%	Once a strf.regtype field exists, the hyperparameters can be updated
%	with nnUpdateHyper
% 
%	E, error due to misestimating the data
%
%
%	Outputs:
%
%	EPRIOR, error due to regularization type
%	E, the total error
%
%Some code taken from Netlab

if isfield(strf,'regtype')
   if strcmp(strf.regtype,'ridge')

        w=nnPak(strf);
	eprior=strf.alpha*(w*w')/2;
   
   elseif strcmp(strf.regtype,'ARD')
	
	w1=strf.w1;
	w1=strf.alpha.*sum(w1.^2,2);
	eprior=sum(w1(:))/2;

   end

   e=e+eprior;
else

eprior=0;

end