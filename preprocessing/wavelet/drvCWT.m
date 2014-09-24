
im=imread('shape.bmp','bmp');
im=double(squeeze(im(:,:,1)));
nlev=5;
wz=cwtDT(im,nlev);
si=cwtCoef2SubIm(wz);

h1=figure;
imagesc(im);
axis image; colormap gray;


[emap,es]=cwtEnergy(si,2,0.25);
%acuMap=cumsum(es,3);
dynRange=max(emap(:))-min(emap(:));
idx=siftlocalmax(emap,dynRange/1000);
[y,x]=ind2sub(size(emap),idx);
keyPos=[x',y'];

figure;
imagesc(emap); axis image; colormap gray; 
axis image; hold on;
plot(x,y,'ko');

ac=cumsum(es,3);
for ii=1:nlev
  emap=es(:,:,ii);
  dynRange=max(emap(:))-min(emap(:));
  idx=siftlocalmax(emap,dynRange/1000);
  [y,x]=ind2sub(size(emap),idx);
  keyPos=[x',y'];

  figure;
  imagesc(emap); axis image; colormap gray; 
  axis image; hold on;
  plot(x,y,'ko');
end  % for ii




[kScale]=cwtKeyScale(emap,keyPos);
figure(h1);
plotCircle(keyPos,kScale);

[pat]=patchExt(im,x',y',2*kScale);

psiz=16;
nLev=log2(psiz);
cwLev=1:nLev;

[normPat]=patchNorm(pat,kScale,kScale,psiz,0);

feature=1  % 1-->16 vector, 2-->60 vector
if feature==1
  rmLev=cwLev(1:end-2)  % remove these
elseif feature==2
  rmLev=[cwLev(1:end-3),cwLev(end)]  % remove these
end
rmLev2=setdiff(cwLev,rmLev)  % Use these 

randPat=magic(psiz);
cwtC=cwtDT(randPat,nLev);
cwtZ=cwtZeroLev(cwtC,rmLev);
wIdx=find(cwtZ.coef~=0);
cwtZ=cwtZeroLev(cwtC,[rmLev,rmLev2]);

clear wdesc wcoef
[wdesc,wcoef]=wdescriptor(normPat,nLev,wIdx);
[descMov,cwtCell]=wdescInv(wdesc,cwtZ,wIdx);

playmovs(1,normPat,descMov)





nPat(ff)=size(normPat{ff},3);













movLen=size(movNorm,3);
for ii=1:movLen
  ii
  im=squeeze(movNorm(:,:,ii));
  wz=cwtDT(im,4);
  si=cwtCoef2SubIm(wz);
  emap=cwtEnergy(si,2,0.2);
  dynRange=max(emap(:))-min(emap(:));
  idx=siftlocalmax(emap,dynRange/1000);
  [y,x]=ind2sub(size(emap),idx);
  keyPos=[x',y'];
  [kScale]=cwtKeyScale(emap,keyPos);

  hf=gcf;
  imagesc(im);
  axis image; colormap gray;
  plotCircle(keyPos,kScale);
  
  drawnow;
  pause(1);
end




