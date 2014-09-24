function [dat,cfStr]=crvPreProc(dat,opt)




% Take curvelet tranform
%--------------------
[dat.eStim,cfStr]=mov2crvlet(dat.eStim,opt.crv);
[dat.vStim,cfStr]=mov2crvlet(dat.vStim,opt.crv);
if ~isempty(dat.cStim)
  [dat.cStim,cfStr]=mov2crvlet(dat.cStim,opt.crv);
end


% Threshold curvelet coefficients
%--------------------
if opt.thresh>0
  dat.eStim=rowThresh(dat.eStim,'ncoef',opt.thresh);
  dat.vStim=rowThresh(dat.vStim,'ncoef',opt.thresh);
  dat.cStim=rowThresh(dat.cStim,'ncoef',opt.thresh);
end



keyboard


