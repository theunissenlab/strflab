function out = df_uppertrivect(in)
%Returns the upper triangular part of the square matrix "in" in a vector.

S = size(in,1);
if S~=size(in,2)
    error(['Error in function uppertrivect: input has size ' num2str(size(in)) ' but it should be a square matrix.']);
end
% 
% ind = zeros(1,(S/2)*(S+1));
% for jj = 1:S
%     ind((jj/2*(jj-1)) + (1:jj)) = (1:jj) +(jj-1)*S;
% end
% out = in(ind)';

% by rows:

out = [];
for jj = 1:S
    out = [out in((jj:S) + (jj-1)*S)];
end
