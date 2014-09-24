function strf=strfUnpak(strf,w)
%function strf=strfUnpak(strf,w)
%
% Separates weights vector into weight and bias matrices and returns a strf 
% data structure identical to the input network, except that the componenet 
% weight matrices have all been set to the corresponding elements of w.  
%
% INPUT:
% [strf] = strf model structure
%    [w] = Parameter vector for STRF, (ex. w=strfPak(strf))
%
% OUTPUT:
% [strf] = strf model structure
%
%
%(Some code modified from NETLAB)

unpakstr=[strf(1).type,'Unpak'];

for ii=1:size(w,1)
  strf(ii)=feval(unpakstr,strf(1),w(ii,:));
end
strf(ii+1:end)=[];
