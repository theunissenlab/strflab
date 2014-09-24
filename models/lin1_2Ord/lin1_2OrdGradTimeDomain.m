function g=lin1_2OrdGradTimeDomain(strf,datIdx)
%function [g] = linGrad(strf, datIdx)
%
% Evaluate gradient of error function for generalized linear model.
%
% Takes a generalized linear model data structure [strf]
% and evaluates the gradient [g] of the error
% function with respect to the network weights. The error function
% corresponds to the choice of output unit activation.  
% 
% 
% INPUT:	
%			[strf] = a linear STRF structure (see linInit for fields)
%			[datIdx] = a set of indices into the global [stim] matrix. 
%
% OUTPUT:
%			   [g] = a 1x(D*L) vector giving the gradient of the error function. 
%					 D=dimensions of stimulus, L=number of delays specified in strf.delays
%
%
%	SEE ALSO
%	linInit, linPak, linUnpak, linFwd, linErr, linGradFourierDomain, linGradTimeDomain
%
%(Some code modified from NETLAB)

global globDat;

samplesize = globDat.nSample;
%  resp_strf = linFwd(strf, 1:globDat.nSample);
[strf,resp_strf] = lin1_2OrdFwd(strf, datIdx);

%  delout = resp_strf - globDat.resp;
if (size(resp_strf,1) == size(globDat.resp(datIdx),2)) & (size(resp_strf,2) == size(globDat.resp(datIdx),1))
    %This means the resps need to be transposed, that's all.
    resp_strf = resp_strf';
end


if strcmp(strf.outputNL, 'huber')
    delout = resp_strf - globDat.resp(datIdx);
    delout = max(-strf.huberD, min(strf.huberD, delout));
elseif strcmp(strf.outputNL, 'logcosh')
    delout = tanh(resp_strf - globDat.resp(datIdx));
else    
    delout = resp_strf - globDat.resp(datIdx);
end

delays1 = strf.delays1;
delays2 = strf.delays2;
numdelays1 = length(delays1);
numdelays2 = length(delays2);
delout(find(isnan(delout))) = 0;

%% assuming nout = 1
delout_array1 = repmat(delout(1)*0, [numdelays1 globDat.nSample]);
delout_array2 = repmat(delout(1)*0, [numdelays2 globDat.nSample]);

for ti=1:numdelays1
	thisIdx = datIdx-delays1(ti);
	validIdx = find(thisIdx>0 & thisIdx<=globDat.nSample);
	delout_array1(ti,thisIdx(validIdx)) = delout(validIdx);
end

for ti=1:numdelays2
	thisIdx = datIdx-delays2(ti);
	validIdx = find(thisIdx>0 & thisIdx<=globDat.nSample);
	delout_array2(ti,thisIdx(validIdx)) = delout(validIdx);
end

gw1 = zeros(size(globDat.stim,2),numdelays1);
gw2 = zeros(size(strf.interactIdx.idxr,1),numdelays2);

%% First order terms

for ii = 1:strf.chunksize:size(globDat.stim,2)
    iiend = min([size(globDat.stim,2) ii+strf.chunksize-1]);
    stim = zeros(samplesize, iiend-ii+1);
    if strf.normalize
        if strf.stds1(ii) ~= 0
            [stim] = norm_std_mean(globDat.stim(strf.covdelays+1:end,ii:iiend), strf.stds1(ii:iiend), strf.means1(ii:iiend));
        else ...

            [stim, strf.stds1(ii:iiend), strf.means1(ii:iiend)] = norm_std_mean(globDat.stim(strf.covdelays+1:end,ii:iiend));
        end
    else
            [stim] = globDat.stim(strf.covdelays+1:end,ii:iiend);
    end
             
    gw1(ii:iiend,:) = (delout_array1*stim)';
end

%% Second order terms

for ii = 1:strf.chunksize:size(strf.interactIdx.idxr,1)
    iiend = min([size(strf.interactIdx.idxr,1) ii+strf.chunksize-1]);
    stim = zeros(samplesize, iiend-ii+1);
    if strf.normalize
        if strf.stds2(ii) ~= 0
            for rr = 1:strf.covdelays + 1
                for cc=1:strf.covdelays + 1
                    idxf = intersect(find(strf.interactIdx.idxrd(ii:iiend)==rr),find(strf.interactIdx.idxcd(ii:iiend)==cc))+ii-1;
                    if ~isempty(idxf)
                        [stim(:,(idxf-ii+1))] = norm_std_mean(globDat.stim((strf.covdelays + 2 - rr):(end-rr+1),strf.interactIdx.idxr(idxf)).*...
                        globDat.stim((strf.covdelays + 2 - cc):(end-cc+1),strf.interactIdx.idxc(idxf)), strf.stds2(idxf), strf.means2(idxf));
                    end
                end
            end
        else ...
            for rr = 1:strf.covdelays + 1
                for cc=1:strf.covdelays + 1
                    idxf = intersect(find(strf.interactIdx.idxrd(ii:iiend)==rr),find(strf.interactIdx.idxcd(ii:iiend)==cc))+ii-1;
                    if ~isempty(idxf)
                        [stim(:,(idxf-ii+1)), strf.stds2(idxf), strf.means2(idxf)] = norm_std_mean(globDat.stim((strf.covdelays + 2 - rr):(end-rr+1),strf.interactIdx.idxr(idxf)).*...
                          globDat.stim((strf.covdelays + 2 - cc):(end-cc+1),strf.interactIdx.idxc(idxf)));
                    end
                end
            end
            

        end
    else ...
        for rr = 1:strf.covdelays + 1
            for cc=1:strf.covdelays + 1
                idxf = intersect(find(strf.interactIdx.idxrd(ii:iiend)==rr),find(strf.interactIdx.idxcd(ii:iiend)==cc))+ii-1;
                if ~isempty(idxf)
                    [stim(:,(idxf-ii+1))] = globDat.stim((strf.covdelays + 2 - rr):(end-rr+1),strf.interactIdx.idxr(idxf)).*...
                      globDat.stim((strf.covdelays + 2 - cc):(end-cc+1),strf.interactIdx.idxc(idxf));
                end
            end
        end
        
    end
    
    gw2(ii:iiend,:) = (delout_array2*stim)';
    
    
    fprintf('.')
end
    
fprintf(' strfGrad step done.\n')
gb1 = nansum(delout);
g = [gw1(:)' gw2(:)'];
if isfield(strf,'b1')
  g = [g gb1];
end
