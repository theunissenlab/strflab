function []=crvCoef(ccoef,mRow,nCol)






% Check Input & Init
%====================
nScale=length(ccoef);
finestLen=length(ccoef{nScale});
if nargin<2


if nargin<2
    coefLen=length(ccoef{nScale});
    if coefLen~=1
        error('crvCoef >< must provide [mRow] & [nCol]');
    else
        coefSize=size(ccoef{nScale}{1});
        mRow=coefSize(1);
        nCol=coefSize(2);
    end
end
coefCell=cell(1,nScale);


[sx,sy,fx,fy,nx,ny]=fdct_wrapping_param(ccoef,mRow,nCol);



for ss=1:nScale
    
    nOri=length(ccoef{ss});
    for oo=1:nOri
        
        [crvX,crvY]=fdct_wrapping_pos2idx(nx,ny,sx,sy,ss,oo,posX,posY);
    end
    
    
    ceofCell{ss}(oo,1:2)=[crvX,crvY];
end





% Extract coef @ scale [scaleLev]
%====================
if scaleLev<=nScale
    scoriCoef=ccoef{scaleLev};
else
    error('crvCoef >< scaleLev must be <= nScale');
end

for ss=1:nSx
    nOx=length(ccx{ss});
    nOy=length(ccy{ss});
    if nOx~=nOy
        error('crvSsim >< # of orientation at each scale must be the same');
    end
    simdex{ss}=zeros(nOx,1);
    
    for oo=1:nOx
        % Process the coefficient at 1 scale and 1 orientation
        cLen=prod(size(ccx{ss}{oo}));
        Cx=reshape(ccx{ss}{oo},cLen,1);
        Cy=reshape(ccy{ss}{oo},cLen,1);
        % Computing Mean
        Mx=mean(Cx);
        My=mean(Cy);
        % Computing Normalized Mean
        nCx=Cx-Mx;
        nCy=Cy-My;
        % Computing 2nd & Product Moments
        Sxy=abs(mean(nCx.*conj(nCy)));
        Sxx=mean(nCx.*conj(nCx));
        Syy=mean(nCy.*conj(nCy));
        % Compute Similarity index at current scale
        %simdex{ss}(oo)=2*(abs(Mx.*conj(My))+k1)*(2*Sxy+k2)./ ...
        %    ((Mx*conj(Mx)+My*conj(My)+k1)*(Sxx+Syy+k2));
        simdex{ss}(oo)=(2*Sxy+k2)./(Sxx+Syy+k2);
    end
end









