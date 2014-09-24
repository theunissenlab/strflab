function [hsp]=crvShowOri(oriCoef,oriRng,hFig)
%function [hsp]=crvShowOri(oriCoef,oriRng,hFig)
% 
% Plot curvelet coefficients at all scale at one orientation.
%
% INPUT:
% [oriCoef] : cell array of curvelet coefficients at all scale with
%             one orientation. oriCoef{s} = coefficients at scale [s].
%    [hFig] : Figure at which the curvelet coeff are plotted.
%
% OUTPUT:
%     [hsp] : subplot handles for the plots.
%
% SEE ALSO: tightsubplot, crvOridx, crvGetOri, crvLet, crvInv
%
% By Michael Wu  --  waftingpetal@yahoo.com (Sep 2006)
%
%====================


% Init Parameters
%====================
if nargin>2
  figure(hFig);
else
  hFig=figure;
end
set(hFig,'Position',[5,5,500,500]);
oriLen=length(oriCoef);


% Plot curvelet coefficients
%====================
fineSize=size(oriCoef{oriLen});
labStr=cell(1,oriLen-1);
if fineSize(1)>fineSize(2)
  for ii=2:oriLen
    hsp(ii)=tightsubplot(1,oriLen-1,ii-1);
    imagesc(abs(oriCoef{ii}));
    axis image; colormap gray; axis off;
    labStr=sprintf(' [ %0.4g , %0.4g )^o',oriRng(ii,1),oriRng(ii,2));
    axlim=axis(gca);
    htx=text(mean(axlim(1:2)),axlim(4),labStr, ...
         'HorizontalAlignment','center','VerticalAlignment','top');
  end
else
  for ii=2:oriLen
    hsp(ii)=tightsubplot(oriLen-1,1,ii-1);
    imagesc(abs(oriCoef{ii}));
    labStr=sprintf(' [ %0.4g , %0.4g )^o',oriRng(ii,1),oriRng(ii,2));
    axis image; colormap gray; axis off;
    axlim=axis(gca);
    htx=text(axlim(2),mean(axlim(3:4)),labStr);
  end
end 

