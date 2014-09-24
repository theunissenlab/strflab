function Cp = LARSrisk(esterr,vars,datIdx)

% function Cp = LARSrisk(esterr,vars,datIdx)
% calculates estimation risk for early stopping in trnLARS
%
% INPUT:
%       esterr = estimation error, calculated in trnLARS
%         vars = number of active variables
%       datIdx = array containing indices for the training set
%
% OUTPUT:
%           Cp = estimation risk
%
% SEE ALSO: trnLARS
% Lucas Pinto, December 2009, lucaspinto@berkeley.edu


global globDat

nonnanidx = find(~isnan(globDat.resp(datIdx)));
Cp = esterr/var(globDat.resp(datIdx(nonnanidx))) - length(nonnanidx) + 2*(vars+1);