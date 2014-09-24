function [astim, stimInfo, params] = preprocSpecVarAGC(stim, filterLen, groupIndex, params)

    grps = unique(groupIndex);
    
    astim = zeros(size(stim));
    P = filterLen;
    
    %% set default parameters
    if ~isfield(params, 'meansub')
        params.meansub = 0;
    end
    if ~isfield(params, 'wmeansub')
        params.wmeansub = 0;
    end
    if ~isfield(params, 'varnorm')
        params.varnorm = 0;
    end
        
    %% set up moving mean filter
    if params.wmeansub
        meanFilter = fliplr(1:P)';
    elseif params.meansub;
        meanFilter = ones(filterLen, 1); %go from tP to t0
    else
        meanFilter = zeros(filterLen, 1); %go from tP to t0        
    end
    
    maxSvar = 5;
    
    %% compute gain-controlled stimulus
    for k = 1:length(grps)
        
        gindx = cv(find(groupIndex == grps(k)));
        tstart = min(gindx);
        for t = gindx        
        
            t0 = max(tstart, t-P+1);
            trng = t0:t;
            shist = stim(trng, :)';            
            mf = flipud(meanFilter(1:length(trng)));
            mfd = max(1, sum(mf));
            
            smean = ((shist * mf) / mfd)';
            svar = 1;
            
            if params.varnorm            
                s2 = bsxfun(@minus, shist', smean).^2;
                svar = mean(s2, 1);                                     
                %maxSvar = max(max(svar), maxSvar);
            end
            
            svgain = 1 - (svar / maxSvar);
            astim(t, :) = (stim(t, :) - smean) .* svgain;
            
        end
                
    end
    maxSvar
        
    stimInfo = struct;
    