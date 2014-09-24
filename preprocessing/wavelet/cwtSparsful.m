function [W]=cwtCoefSparsful(W)
%function [W]=cwtCoefSparsful(W)
%
% Convert between sparse and full CWT coefficient in the .coef
% field of CWT coefficient structure. If [W] is a CWT structure,
% it will do the conversion. If it is a full matrix, it will 
% convert it to sparse. If it is a sparse matrix, it will convert
% it to full. If [W] is a cell of CWT structure, then it will 
% call it self again with each elements of the cell as input.
%
% INPUT:
% [W] = CWT structure, with [.coef] field and [.siz] field. Or it
%       can be a cell array of CWT structures.
% OUTPUT:
% [W] = Same as input, except the [.coef] field is converted to 
%       full or sparse depending on what the input was.
%
% SEE ALSO: wdescriptor, wdescHist, wdescInv, cwtDT, cwtInv
%
% By Michael Wu  --  waftingpetal@yahoo.com (Apr 2007)
%
% ====================


if isstruct(W)
  if issparse(W.coef)
    W.coef=full(W.coef);
  else
    W.coef=sparse(W.coef);
  end  % issparse
elseif iscell(W)
  nCC=length(W);
  if issparse(W{1}.coef)
    for ii=1:nCC
      W{ii}.coef=full(W{ii}.coef);
    end  % for ii
  else
    for ii=1:nCC
      W{ii}.coef=sparse(W{ii}.coef);
    end  % for ii
  end  % issparse
end  % if isstruct


%  for ii=1:nCC
%    W{ii}=cwtSparsful(W{ii});
%  end  % for ii




