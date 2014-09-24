function [strf,resp_strf_dat,im_dat,fm_dat,hl_dat,a_dat]=cnnFwd(strf,datIdx)

%NNFWD	Forward propagation through 2-layer network.
%
%	Description
%	RESP_STRF=NNFWD(STRF, STIM) takes a neural network model data 
%	structure STRF together with a stimulus matrix STIM and forward 
%	propagates the inputs through the network to generate a vector 
%	RESP_STRF of output vectors.
%
%	[RESP_STRF, Z, A] = NNFWD(STRF, STIM) also returns a matrix STIM 
%	giving the summed inputs to each hidden unit in A and out of
%	each hidden unit in Z unit, where each row corresponds to one 
%	pattern.
%
%
%	Usage:
%
%	[RESP_STRF, A, Z]=NNFWD(STRF, STIM)
%
%
%	Inputs:
%
%	STRF,	A STRF structure
%	STIM,	An NXM stimulus matrix (or vector)
%
%
%	Outputs:
%
%	RESP_STRF, Predicted output elicited from STRF structure by
%		   STIM.
%	A, 	Summed inputs to each hidden unit
%	Z,	Hidden unit activations
%
%
%	See also
%	NNINIT, NNPAK, NNUNPAK, NN1ERR, NNGRAD
%
%
%(Some code modified from NETLAB)

global globDat;

samplesize = globDat.nSample;

if strf.internal.compFwd == 0 & samplesize == length(strf.internal.prevResp) & strf.internal.dataHash == globDat.dataHash
	resp_strf_dat = strf.internal.prevResp(datIdx);
	a_dat = strf.internal.prevLinResp(datIdx);
        fm_dat = strf.internal.prevHidResp(datIdx);
	return
end

% z = zeros(samplesize, strf.nHidden);
fm = zeros((strf.frameSize(1)-strf.filtDim(1)+1),(strf.frameSize(2)-strf.filtDim(2)+1), size(strf.filts,3),samplesize);
for fi = 1:size(strf.filts,3)
    for ti=1:length(strf.delays)
        
        if (strf.constrain == 1 | strf.constrain == 2) & ti == 1
            zt(:,:,:,fi,ti) = convn(globDat.stim, squeeze(strf.filts(:,:,fi)), 'valid');
            tii = 1;
        elseif (strf.constrain == 0 | strf.constrain == 3)
            zt(:,:,:,fi,ti) = convn(globDat.stim, squeeze(strf.filts(:,:,fi,ti)), 'valid');
            tii = ti;
        else
            tii = 1;
        end
        thisshift = strf.delays(ti);
        if thisshift>=0            
            fm(:,:,fi,thisshift+1:end) = squeeze(fm(:,:,fi,thisshift+1:end)) + strf.fmw(fi,ti)*zt(:,:,1:end-thisshift,fi,tii);
            
        else
            offset = mod(thisshift, samplesize);
            fm(:,:,fi,1:offset) = squeeze(fm(:,:,fi,1:offset)) + strf.fmw(fi,ti)*zt(:,:,-thisshift+1:end,fi,tii);
        end
    end
end

if ismember(strf.constrain, [1 3])
    for ii = 1:size(z,3)
        fm(:,:,ii,:) = tanh(fm(:,:,ii,:) + strf.fmb(ii));
    end
else
    fm = tanh(sum(fm,3) + strf.fmb);
end
         
% 
% for ti=1:length(strf.delays)
%   zt = globDat.stim * strf.w1(:,:,ti);
% 
%   thisshift = strf.delays(ti);
%   if thisshift>=0
%     z(thisshift+1:end,:) = z(thisshift+1:end,:) + zt(1:end-thisshift,:);
%   else
%     offset = mod(thisshift, samplesize);
%     z(1:offset,:) = z(1:offset,:) + zt(-thisshift+1:end,:);
%   end
% 
% end


% z = z + repmat(strf.b1,[size(z,1) 1]);

% mask for nonvalid frames
nanmask = mod(strf.delays, size(globDat.stim,1)+1);
nanmask = nanmask(find(nanmask)); % no mask for delay 0
fm(:,:,:,nanmask) = NaN;
% z=tanh(z);
szz = size(fm);
fm = reshape(fm, [prod(szz(1:3)) szz(4)])';

if ~isempty(strf.hlw)
    hl = fm * strf.hlw;
    hl = hl + repmat(strf.hlb, [size(a,1) 1]);
    hl = tanh(hl);
else
    hl = fm; 
end

% clear z

a = hl*strf.outw + strf.outb;

switch strf.outputNL

  case 'linear'     % Linear outputs
    resp_strf = a;

  case 'logistic'   % Logistic outputs
    % Prevent overflow and underflow: use same bounds as glmerr
    % Ensure that log(1-y) is computable: need exp(a) > eps
    maxcut = -log(eps);
    % Ensure that log(y) is computable
    mincut = -log(1/realmin - 1);
    a = min(a, maxcut);
    a = max(a, mincut);
    resp_strf = 1./(1 + exp(-a));

  case 'softmax'        % Softmax outputs
    nout = size(a,2);
    % Prevent overflow and underflow: use same bounds as glmerr
    % Ensure that sum(exp(a), 2) does not overflow
    maxcut = log(realmax) - log(nout);
    % Ensure that exp(a) > 0
    mincut = log(realmin);
    a = min(a, maxcut);
    a = max(a, mincut);
    temp = exp(a);
    resp_strf = temp./(sum(temp, 2)*ones(1,nout));
    % Ensure that log(y) is computable
    resp_strf(resp_stf<realmin) = realmin;

  case 'exponential'
    resp_strf=exp(a);
  otherwise
    error(['Unknown activation function ', strf.outputNL]);
end


resp_strf_dat = resp_strf(datIdx);
a_dat = a(datIdx);
fm_dat = fm(datIdx,:);
im_dat = zt(:,:,datIdx,:,:);
if isempty(strf.hlw)
    hl_dat = [];
else
    hl_dat = hl(datIdx,:);
end

%strf.internal.compFwd = 0;
%strf.internal.prevResp = resp_strf;
%strf.internal.prevLinResp = a;
%strf.internal.prevHidResp = z;
%strf.internal.dataHash = globDat.dataHash;