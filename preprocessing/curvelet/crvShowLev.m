function [hsp]=crvShowLev(levCoef,hFig)
%function [hsp]=crvShowLev(levCoef,hFig)
% 
% Plot all curvelet coefficients at all orientation at one scale.
%
% INPUT:
% [levCoef] : cell array of curvelet coefficients of all orientations at
%             one scale. scaleCoef{w} = coefficient at orientation [w].
%    [hFig] : Figure at which the curvelet coeff are plotted.
%
% OUTPUT:
%     [hsp] : subplot handles for the plots.
%
% SEE ALSO: tightsubplot, crvOridx, crvGetLev, crvLet, crvInv
%
% By Michael Wu  --  waftingpetal@yahoo.com (Sep 2006)
%
%====================


% Init Parameters
%====================
if nargin>=2 & ishandle(hFig);
  figure(hFig);
else
  hFig=figure;
end
set(hFig,'Position',[5,5,500,500]);
nOri=length(levCoef);


% Plot Polar Grid
%====================
nsq=nOri/4;
sideLen=2+round(nsq);

tightsubplot(1,1,1);
axis([-1,1,-1,1]); axis off;
polGrid=(3*pi/4)-linspace(0,2*pi,nOri+1);
zeroGrid=zeros(size(polGrid));
oneGrid=zeroGrid+1;

[posX,posY]=pol2cart(polGrid,oneGrid*2);
line([zeroGrid;posX],[zeroGrid;posY],'color','k');

pGridSiz=polGrid(2)-polGrid(1);
txtGrid=polGrid(1:nOri)+pGridSiz/2;
[posX,posY]=pol2cart(txtGrid,oneGrid(1:nOri)*0.6);
for ii=1:nOri
  htx=text(posX(ii),posY(ii),num2str(ii),'HorizontalAlignment','center');
end


% Plot curvelet coefficients
%====================
hsp=zeros(nOri,1);
if nOri==1
  hsp=imagesc(abs(levCoef{1}));
  axis off; axis image; colormap gray;
else
  for ii=1:nsq
    hsp(ii)=tightsubplot(sideLen,sideLen,ii+1);
    imagesc(real(levCoef{ii}));
    axis off; axis image; colormap gray;
        
    hsp(ii+nsq)=tightsubplot(sideLen,sideLen,(ii+1)*sideLen);
    imagesc(real(levCoef{ii+nsq}));
    axis off; axis image; colormap gray;

    hsp(ii+2*nsq)=tightsubplot(sideLen,sideLen,(sideLen^2)-ii);
    imagesc(real(levCoef{ii+2*nsq}));
    axis off; axis image; colormap gray;

    hsp(ii+3*nsq)=tightsubplot(sideLen,sideLen,(nsq+1-ii)*sideLen+1);
    imagesc(real(levCoef{ii+3*nsq}));
    axis off; axis image; colormap gray;
  end
end


