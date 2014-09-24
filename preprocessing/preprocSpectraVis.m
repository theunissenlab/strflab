function preprocSpectraVis(net);
%function preprocSpectraVis(net);
%
% A visualizer of glm net preprocessed by preprocSpectra
% 
% INPUT:
% [net] = strf structure to be visualized
%

params = net.params;

fsize = params.fSize;
w = net.w1;
delays = net.delays;

k = reshape(w, [fsize length(delays)]);

maxk = max(abs(k(:)));

for pi=1:fsize(4)
	figure(pi); clf;
	for ti=1:fsize(3)
		for di=1:length(delays)
			subplot(length(delays), fsize(3), ti+(di-1)*fsize(3));
			thisim = squeeze(k(:,:,ti,pi,di));
			imagesc(thisim, [-maxk maxk]);
			axistickoff();
			if ti==1, ylabel(sprintf('d: %d', delays(di))), end
		end
	end
end

return;


function axistickoff()

set(gca, 'xtick', []);
set(gca, 'ytick', []);

return;
