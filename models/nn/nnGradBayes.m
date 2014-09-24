function [g, gprior]=nnGradBayes(strf, g)
%NNERR	Evaluate error function for 2-layer network.
%
%	Description
%	[E, EPRIOR]= NNGRADBAYES(STRF, G) takes a network data structure STRF
%	and a gradient vector G, and adds the gradient due to regularization.
%	
%
%	Inputs:
%	
%	STRF, a neural network strf structure
%		nnGradBayes uses the field strf.regtype to determine how
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
%	G, gradient due to misestimating the data
%
%
%	Outputs:
%
%	GPRIOR, gradient due to regularization type
%	G, the total gradient
%
%Some code taken from Netlab
if isfield(strf,'regtype')
   if strcmp(strf.regtype,'ridge')

        w=nnPak(strf);
	gprior=strf.alpha*w;
   
   elseif strcmp(strf.regtype,'ARD')
	
	strf.w1=repmat(strf.alpha,[1 strf.nhidden 1]).*strf.w1;
        strf.b1=strf.b1*0;
	strf.w2=strf.w2*0;
	strf.b2=strf.b2*0;
	gprior=nnPak(strf);

   end

   g=g+gprior;
else

   gprior=0*g;
end