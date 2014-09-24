function y = randsample(pop, k, varargin)
%RANDSAMPLE Random sample, with or without replacement.
%   Y = RANDSAMPLE(N,K) returns Y as a vector of values sampled uniformly
%   at random, without replacement, from the integers 1:N.
%
%   Y = RANDSAMPLE(POPULATION,K) returns K values sampled uniformly at
%   random, without replacement, from the values in the vector POPULATION.
%
%   Y = RANDSAMPLE(...,REPLACE) returns a sample taken with replacement if
%   REPLACE is true, or without replacement if REPLACE is false (the default).

if (length(varargin) > 0 && varargin{1})
% Replace is TRUE
	if (length(pop) == 1)
		y = ceil(rand([k, 1]).*pop);
	else
		y = pop(ceil(rand([k,1]).*length(pop)));
	end
else
	if (length(pop) == 1)
		y = randperm(pop);
	else
		y = pop(randperm(length(pop)));
	end
	y = y(1:min(k, length(y)));
end
