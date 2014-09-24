function g=nnGradTimeDomain(strf,datIdx)

%NNGRAD Evaluate gradient of error function for generalized linear model.
%
%	Description
%	G = NNGRAD(STRF, STIM, RESP) takes a neural networ generalized data
%	structure STRF  together with a matrix STIM of input vectors and a 
%	matrix RESP of responses, and evaluates the gradient G of the error
%	function with respect to the network weights. The error function
%	corresponds to the choice of output unit activation.
%
%
%	Usage:
%	
%	G = NNGRAD(STRF, STIM, RESP)
%
%
%	Inputs:
%		
%	STRF,	A STRF structure
%	STIM,	An nXm stimulus matrix (or vector)
%	RESP,	An nXtrials response matrix (or vector)
%
%
%	Outputs:
%
%	G, Value of the gradient of the error function for a given STRF,
%	   STIM and RESP
%
%
%
%	See also
%	NNINIT, NNPAK, NNUNPAK, NNFWD, NNERR
%
%(Some code modified from NETLAB)

global globDat;


%  resp_strf = nnFwd(strf, 1:globDat.nSample);
[strf,resp_strf,z,a] = nnFwd(strf, datIdx);

%  delout = resp_strf - globDat.resp;
if (size(resp_strf,1) == size(globDat.resp(datIdx),2)) & (size(resp_strf,2) == size(globDat.resp(datIdx),1))
    %This means the resps need to be transposed, that's all.
    resp_strf = resp_strf';
end


delout = resp_strf - globDat.resp(datIdx);

delays = strf.delays;
numdelays = length(delays);
delout(find(isnan(delout))) = 0;
z(find(isnan(z))) = 0;


gw2 = z'*delout;
gb2 = sum(delout, 1);

delhid = delout*strf.w2';
delhid = delhid.*(1.0 - z.*z);

gb1 = sum(delhid, 1);

g=[];

for ii=1:strf.nHidden

  %% assuming nout = 1
     delout_array = zeros(numdelays,globDat.nSample);

  for ti=1:numdelays
	thisIdx = datIdx-delays(ti);
	validIdx = find(thisIdx>0 & thisIdx<=globDat.nSample);
	delout_array(ti,thisIdx(validIdx)) = delhid(validIdx,ii);
  end

  gw1 = (delout_array*globDat.stim)';


  g=[g gw1(:)'];

end

g = [g gb1];

g= [g gw2'];

g= [g gb2];



