function out = displayAuditoryStrf(strf);
filter = squeeze(strf.w1);
out1 = simagesc(filter);
if nargout > 0
    out = out1;
end