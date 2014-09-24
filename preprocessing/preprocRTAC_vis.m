function [filtRTA, filtRTC1, filtRTC2] = preprocRTAC_vis(strf)

%function [out, strf] = preprocRTAC_Vis(strf)
%
%  A visualizer of STRF preprocessed by preprocRTAC
%
% INPUT:
%          [strf] = strf model that was preprocessed by preprocWavelets3d
%
% OUTPUT:
%       [filtRTA] = a matrix of the RTA filter
%      [filtRTC1] = a matrix of the 1st Eigenvector of the RTC (1st RTC filter)
%      [filtRTC1] = a matrix of the 2nd Eigenvector of the RTC (2nd RTC filter)
%


if strf.params.covtime
    delays = (strf.params.covdelays + 1);
else
    delays = length(strf.delays);
end
lines = 0;
covstart = 1;
if strf.params.RTAC(1)
    RTA = strf.w1(1:(strf.params.covdelays + 1)*strf.params.covsz,:);
    lines = 1;
    covstart = (strf.params.covdelays + 1)*strf.params.covsz+1;
else
    RTA = zeros([strf.params.covsz delays]);
end

if strf.params.RTAC(2)
    lines = lines + 2;
    RTCsub = strf.w1(covstart:end,:);

    if strf.params.covtime
        idx2 =sub2ind([strf.params.covsz*(strf.params.covdelays + 1) strf.params.covsz*(strf.params.covdelays + 1)], strf.params.idxr, strf.params.idxc);
        RTC = squeeze(zeros([strf.params.covsz*(strf.params.covdelays + 1) strf.params.covsz*(strf.params.covdelays + 1)]));
        RTC(idx2) = RTCsub;
        [u,s,v] = svd(RTC);
        delays = (strf.params.covdelays + 1);
    else
        
        idx2 =sub2ind([strf.params.covsz strf.params.covsz], strf.params.idxr, strf.params.idxc);

        for ii = 1:delays
            RTC = squeeze(zeros([strf.params.covsz strf.params.covsz]));
        
            RTC(idx2) = RTCsub(:,ii);

            RTC = RTC + tril(RTC', -1);

            [u(:,:,ii),s(:,:,ii),v(:,:,ii)] = svd(RTC);

        end
    end
else
    u = zeros([strf.params.covsz 2 delays]);
end

filtRTA = reshape(RTA, [sqrt(strf.params.covsz) sqrt(strf.params.covsz) delays]);
filtRTC1 = reshape(u(:,1,:), [sqrt(strf.params.covsz) sqrt(strf.params.covsz) delays]);
filtRTC2 = reshape(u(:,2,:), [sqrt(strf.params.covsz) sqrt(strf.params.covsz) delays]);

%  Plot results -----------
if strf.params.RTAC(2)
for yy = 1:size(s,3)
figure;
plot(diag(s(:,:,yy)), 'o'); % examine eigenvalues
title('eigenvalues');
end
end

RTAmax = max([abs(filtRTA(:)); 1e-6]);
RTC1max = max([abs(filtRTC1(:)); 1e-6]);
RTC2max = max([abs(filtRTC2(:)); 1e-6]);
f2 = figure;
set(f2, 'Position', [150 150 1000 500]);
colormap gray;

for iii = 1:delays
    if strf.params.RTAC(1)
        subplot(lines,delays,iii); imagesc(filtRTA(:,:,iii), [-RTAmax RTAmax]); axis image; axis off; title(['RTA: ' num2str(-iii+1)]);
        nextline = delays;
    else
        nextline = 0;
    end
    if strf.params.RTAC(2)
        subplot(lines,delays,iii+nextline); imagesc(filtRTC1(:,:,iii), [-RTC1max RTC1max]); axis image; axis off; title(['RTC1: ' num2str(-iii+1)]);
        subplot(lines,delays,iii+nextline+delays); imagesc(filtRTC2(:,:,iii), [-RTC2max RTC2max]); axis image; axis off; title(['RTC2: ' num2str(-iii+1)]);
    end
end
