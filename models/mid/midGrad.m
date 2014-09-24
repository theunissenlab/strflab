function [strf,g,savedStimFft]=midGrad(strf,datIdx,savedStimFft)
%function [strf, g,savedStimFft] = midGrad(strf, datIdx, savedStimFft)
%
% Evaluate gradient of error function for generalized linear model.
%
% Takes a generalized linear model data structure [strf]
% and evaluates the gradient [g] of the error
% function with respect to the network weights. The error function
% corresponds to the choice of output unit activation.  If [savedStimFft]
% is included as an input, the gradient will be evaluated using a fourier domain 
% technique.  Alternatively, this method may be chosen by setting strf.fourierdoman = 1.
% 
% 
%
% INPUT:	
%			[strf] = a linear STRF structure (see linInit for fields)
%			[datIdx] = a set of indices into the global [stim] matrix. 
%   [savedStimFft] = (optional) the Fourier transform of the stimulus, which is
%   				 sometimes saved and which makes Fourier domain calculations of the grad
%   				 much faster.
%
% OUTPUT:
%			[strf] = unmodified STRF structure.
%			   [g] = a 1x(D*L) vector giving the gradient of the error function. 
%					 D=dimensions of stimulus, L=number of delays specified in strf.delays
%   [savedStimFft] = (optional) if the fourier-based gradient has been selected, will contain
%					 the Fourier transform of the stimulus.
%
%
%	SEE ALSO
%	linInit, linPak, linUnpak, linFwd, linErr, linGradFourierDomain, linGradTimeDomain
%
%(Some code modified from NETLAB)

global globDat;

% Check arguments for consistency
% errstring = consist(strf, 'mid', x, t);
% if ~isempty(errstring);
%   error(errstring);
% end

% y = midfwd(strf, x);
[strf,resp_strf] = midFwd(strf, datIdx);

delays = strf.delays;
numdelays = length(delays);

x=[-strf.nbin*strf.bin:strf.bin:strf.nbin*strf.bin];
y=resp_strf;
s=globDat.stim(datIdx,:);
r=globDat.resp(datIdx);

px=histc(y,x);
px(:)=px(:)/sum(px);

z=find(r>0);
pxr=histc(y(z),x);
pxr(:)=pxr(:)/sum(pxr);

sb=zeros(size(s,2),size(x,2), numdelays);
sbr=zeros(size(s,2),size(x,2), numdelays);

for ii = max([find(isnan(y))+1; 1]):size(s,1)
% for ii=1:size(s,1)
  bx=find(histc(y(ii),x)>0);
	for di = 1:numdelays
		if ii-di+1 >= 1
  			sb(:,bx,di)=sb(:,bx,di)+s(ii-di+1,:)';
  			if(r(ii)>0)
    			sbr(:,bx,di)=sbr(:,bx,di)+s(ii-di+1,:)';
  			end;
		end
	end
end;

del=zeros(size(s,2), numdelays);
a=zeros(size(s,2), numdelays);
z=find(pxr>0);

for di2 = 1:numdelays
	for ii=z(1):z(end)
  
  		%ii,px(ii+1),px(ii-1)
  		if((px(ii+1)>0)&(px(ii-1)>0))
    	%keyboard;
    	a(:,di2)=pxr(ii)*(sbr(:,ii, di2)'-sb(:,ii, di2)')*(pxr(ii+1)/px(ii+1)-pxr(ii-1)/px(ii-1))/strf.bin;
    	%keyboard;
    	del(:,di2)=del(:,di2)+a(:,di2);
    	%del
  		end;
	end;
end

gw1=-del;

% gdata=g;
% delays = strf.delays;
% numdelays = length(delays);
% delout(find(isnan(delout))) = 0;
% 
% %% assuming nout = 1
% delout_array = repmat(delout(1)*0, [numdelays globDat.nSample]);
% 
% for ti=1:numdelays
% 	thisIdx = datIdx-delays(ti);
% 	validIdx = find(thisIdx>0 & thisIdx<=globDat.nSample);
% 	delout_array(ti,thisIdx(validIdx)) = delout(validIdx);
% end

% gw1 = (delout_array*globDat.stim)';
% gb1 = nansum(delout);
g = gw1(:)';
% if isfield(strf,'b1')
%   g = [g gb1];
% end
g = [g 0];

%[g, gdata, gprior] = gbayes(strf, gdata);
