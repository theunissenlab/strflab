function [preds, maxidx] = optimizeLARSstep(strfTrained, options, cvResult)

strf = strfTrained(1);

for ii = size(options,2)
    steps(ii) = size(options(ii).diagnostics.betas,2);
end



for ii = 1:size(options,2)
    ii
	for iii=1:max(steps)
	    
	    if iii <= size(options(ii).diagnostics.betas,2)
            strf = strfUnpak(strf,[options(ii).diagnostics.betas{iii} strfTrained(ii).b1]);
            [strf,predResp]=strfFwd(strf,cvResult(ii).testIdx);
        
            nonnanidx = find(~isnan(predResp));
            
            if ~isempty(nonnanidx)
                preds(ii, iii) = corr(predResp(nonnanidx), cvResult(ii).true(nonnanidx));
            else
                preds(ii, iii) = 0;
            end
            
        else
            preds(ii, iii) = 0;
        end
    end
end

plot(preds');
maxpreds = max(preds');

for vv = 1:size(preds,1)
    maxidx(vv) = min(find(preds(vv,:) == maxpreds(vv)));
    hold on; plot(maxidx(vv), preds(vv,maxidx(vv)), 'x', 'Color', 'k');
    text(double(maxidx(vv)), double(preds(vv,maxidx(vv))), ['(' num2str(maxidx(vv)) ', ' num2str(preds(vv,maxidx(vv))) ')']);
end







