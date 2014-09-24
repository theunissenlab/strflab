%% Convolve a stimulus with a STRF
%
%   Input:
%       allstim:
%
%       delays:
%
%       strf:
%
%       groupIndex:
%
%   Output:
%
%       modelResponse:
%
%
function modelResponse = conv_strf(allstim, delays, strf, groupIndex)

    ugi = unique(groupIndex);
    nDatasets = length(ugi);
    timeLen = size(allstim, 1);
    a = zeros(timeLen, 1);
    
    for k = 1:nDatasets
        
        gi = ugi(k);
        rng = find(groupIndex == gi);
        soff = rng(1) - 1;
        stim = allstim(rng, :);
        for ti = 1:length(delays)
            at = stim * strf(:, ti);
            eindx = soff + length(at);
            
            thisshift = delays(ti);
            if thisshift >= 0                
                a(soff+thisshift+1:eindx) = a(soff+thisshift+1:eindx) + at(1:end-thisshift);
            else
                offset = eindx + thisshift;
                a(soff+1:offset) = a(soff+1:offset) + at(-thisshift+1:end);
            end
        end
    end
    
    a(isnan(a)) = 0;
    modelResponse = a';
