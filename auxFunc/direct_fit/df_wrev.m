function y = df_wrev(x)
%WREV Flip vector.
%   Y = WREV(X) reverses the vector X.
%
%   See also FLIPLR, FLIPUD.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 01-May-96.
%   Last Revision: 01-Jun-1998.
%   Copyright 1995-2000 The MathWorks, Inc.
% $Revision: 1.8 $

% Check arguments.
if df_errargn(mfilename,nargin,[1],nargout,[0:1]), error('*'), end

y = x(end:-1:1);
