function [s_stds, s_means] = norm_std_mean(s_stds, s_means);

global s

divnum = 20;
snum = size(s,2);
chunksize = ceil(snum/divnum);

if ~exist('s_means', 'var')
	%% calculate mean
	s_means = mean(s, 1);
end

if ~exist('s_stds', 'var')
	%% calculate std
	if prod(size(s)) < 100000 | size(s,2) < 20
		s_stds = std(s, 0, 1);
	else % to avoid out of memory, calculate stds by chunks
		st = 1;
		for ii=1:divnum
			ed = min([snum st+chunksize-1]);
			s_stds(st:ed) = std(s(:,st:ed), 0, 1);
			st = st + chunksize;
		end
	end
end


%% offset for mean
if prod(size(s)) < 100000 | size(s,2) < 20
	meanmat = repmat(s_means, [size(s,1) 1]);
	s = s - meanmat;
else
	st = 1;
	for ii=1:divnum
		ed = min([snum st+chunksize-1]);

		meanmat = repmat(s_means(st:ed), [size(s,1) 1]);
		s(:,st:ed) = s(:,st:ed) - meanmat;

		st = st + chunksize;
	end
end


s_stds(find(s_stds==0)) = 1; % To avoid zero-division

%% normalize stds
if prod(size(s)) < 100000 | size(s,2) < 20
	stdmat = repmat(s_stds, [size(s,1) 1]);
	s = s./stdmat;
else
	st = 1;
	for ii=1:divnum
		ed = min([snum st+chunksize-1]);

		stdmat = repmat(s_stds(st:ed), [size(s,1) 1]);
		s(:,st:ed) = s(:,st:ed)./stdmat;

		st = st + chunksize;
	end
end
