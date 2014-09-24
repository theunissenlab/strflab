function strfData(stim,resp,groupIdx)
%function strfData(stim,resp,groupIdx)
%
% Takes [stim] and [resp] from the preprocessing routines, and set them
% as global data.
%
% INPUT:
%     [stim] = NxD matrix of N samples and D dimensions.
%     [resp] = Nx1 column vector of N samples.
% [groupIdx] = Nx1 vector of integers, specifying which group the [stim] samples
%              belongs. If it is empty, each sample is its own group.
% OUTPUT: 
% No output, but sets the following global variable. You must have a line
%   global globDat;
% in your code in order to use the global variable [globDat].
%
%  [globDat] = global variable containing following fields:
%      .stim = same as input [stim]
%      .resp = same as input [resp]
%   .nSample = # of sampels = N, the common dimension of [stim] and [resp].
%
% SEE ALSO:
%
% By Michael Wu  --  waftingpetal@yahoo.com (Aug 2007)
% Edited by Michael Oliver to remove slow STRFPAK hashes(May 2009)
%
% ====================



% Init & get input size
% --------------------
stimSiz=size(stim);
respLen=length(resp);
nSample=intersect(stimSiz,respLen);
if length(nSample)~=1
  warning('strfData >< [stim] and [resp] must have same # of rows.');
end
if nargin>2
  idxLen=length(groupIdx);
  if idxLen~=nSample
    warning('strfData >< [groupIdx] must have same length as [resp].');
  end
else
  groupIdx=[];
end


% Init & get input size
% --------------------
global globDat;
globDat.stim=stim;
globDat.resp=resp;
globDat.nSample=nSample;
globDat.groupIdx=groupIdx;


% Compute hash for the data set 
% --------------------
hresp = resp;
if iscell(resp)
    hresp = [];
    %resp is a cell array of spike time vectors
    for k = 1:length(resp)
        st = resp{k};
        hresp = [hresp st(:)'];
    end
end

respHash = 100*abs(nanmean(double(hresp(:))) + nanmean(double(hresp(1:11:end))));
stimHash = 100*abs(nanmean(double(stim(1:109:end))));
magdif = log10((respHash+0.00012)/(stimHash+0.00011));
dataHash = respHash + stimHash * 10^magdif;
globDat.dataHash = dataHash;




