
% Get V4 neuron with Natural image sequence in pixel/psth format
%--------------------
close all;
clear all;
[dbDat,cellList,cellDBidx,nCells]=getDBparm('v2',2,1,0);
nFrame=500  % 0=all frames
nDim=48;


% Get neuron logistics 
%--------------------
stimSizes=cell2mat({dbDat.stimwindowsize});
stimCRFs=cell2mat({dbDat.stimfilecrf});
nTrial=cell2mat({dbDat.repcount});
stimRate=cell2mat({dbDat.stimspeedid});
nspikes=cell2mat({dbDat.spikes});
resplen=cell2mat({dbDat.resplen});
meanRate=nspikes./resplen;


% Pick a few good neurons
%--------------------
testId=intersect(find(nTrial==1), ...
	intersect(find(meanRate>1.3), ...
	intersect(find(stimSizes>nDim),find(stimCRFs>1))));
if 0
testId=intersect(find(meanRate>1), ...
	intersect(find(stimSizes>nDim),find(stimCRFs>1)));
end  %if 0;


% found cells with mean firing rate >1, stim >64 pixel, and >1CRF
%--------------------
if 0
stimSizes(testId)
stimCRFs(testId)
nTrial(testId)
meanRate(testId)


% Pick the data set to analyze
%--------------------
testDat=testId(3);
cellId=dbDat(testDat).cellid
end  % if 0
goodCell={'e0017','e0024','e0047','e0100'};
cellId=goodCell{4}
%cellNo=strmatch('e0022',cellList)
cellNo=strmatch(cellId,cellList);

cellDBidx{cellNo}
stimSizes(cellNo)
stimCRFs(cellNo)
nTrial(cellNo)
stimRate(cellNo)
nspikes(cellNo)
resplen(cellNo)
meanRate(cellNo)


% Get stim covering 1CRF 64x64 pixel, response in psth
%--------------------
parm=dbDat(cellDBidx{cellNo});
%parm=dbDat(testDat);
parm=setStimParm(parm,nDim,4.5);
parm=setRespParm(parm,0);

tic;
nlag=5
[dat,parm]=cellStimResp(parm,0.1,0.1,nFrame);
dat=preProcData(dat);
%[estDat.stim,estDat.resp]=tapDelNan(estDat.stim,estDat.resp,nlag);    
%[estDat.cvStim,estDat.cvResp]=tapDelNan(estDat.cvStim,estDat.cvResp,nlag);
%[valDat.stim,valDat.resp]=tapDelNan(valDat.stim,valDat.resp,nlag); 
if nFrame==0
  nFrame=size(dat.eStim,3);
end
t.getData=toc; time2str(toc);

%movVec=estDat.stim';
%resp=estDat.resp;



% clear double frames
%--------------------
frameInc=60/parm.stimspeedid;
if frameInc>1
  rmIdx=1:frameInc:nFrame;
  movVec(:,rmIdx)=[];
  nFrame=size(movVec,2);
end


% Normalize movie frame to [0,1] for SIFT computation
%--------------------
movMin=min(movVec);
movVec=movVec-repmat(movMin,nDim^2,1);
movMax=max(movVec);
movNorm=reshape(movVec./repmat(movMax,nDim^2,1),nDim,nDim,nFrame);
clear movVec;  % Don't need it anymore


% SIFT key point extraction
%--------------------
addpath('/auto/k1/wafting/mat/pack/sift/');
%addpath('/auto/k2/share/matlab/ivanov/sift/');
%addpath('/auto/k2/share/matlab/ivanov/siftutils/')
nSubLev=3;
keyPt=cell(1,nFrame);

clear keyPt;
nKey=zeros(1,nFrame);
keyPt=cell(1,nFrame);
tic
for ff=1:nFrame  % 1:nFrame
  [keyPt{ff}]=swift(movNorm(:,:,ff),'Threshold',0.003, ...
  'NumLevels',nSubLev,'FirstOctave',-1);
  nKey(ff)=size(keyPt{ff},2);
  
  if 0  
  figure;
  imagesc(movNorm(:,:,ff));
  colormap gray;axis image;
  plotsiftframe(keyPt{ff});  % See Figure
  drawnow;
  end  % if 0
end
hist(nKey,floor(sqrt(sum(nKey))));
nzKey=find(nKey>0);
t.sift=toc; time2str(toc);
%figure;plotss(gss{ff});colormap gray
%figure;plotss(dogss);colormap gray


% Prepare for extract keyPt Patches & Normalized patches
%--------------------
addpath('/auto/k1/wafting/mat/pack/cwt/');
%addpath /auto/k1/wafting/mat/pack/pyrTools/;
psiz=16;
nLev=log2(psiz);
randPat=magic(psiz);
cwtC=cwtDT(randPat,nLev);
%cwtSubImg(cwtC)


% Pick a feature descriptor
%--------------------
feature=2  % 1-->16 vector, 2-->60 vector
if feature==1
  rmLev=[1 2 3];
  rmLev2=[4 5];
elseif feature==2
  rmLev=[1 2 5];
  rmLev2=[3 4];
end
cwtZ=cwtZeroLev(cwtC,rmLev);
wIdx=find(cwtZ.coef~=0);
cwtZ=cwtZeroLev(cwtC,[rmLev,rmLev2]);


% Extraction and Normalization
%--------------------
nPat=zeros(1,nFrame);
doFrame=nFrame;
normPat=cell(1,nFrame);
tic
for ff=nzKey
  [pat]=patchExt(movNorm(:,:,ff),keyPt{ff}(1,:),keyPt{ff}(2,:),5*keyPt{ff}(3,:));
  [normPat{ff}]=patchNorm(pat,keyPt{ff}(3,:),keyPt{ff}(4,:),psiz,0);
  nPat(ff)=size(normPat{ff},3);
end
nzPat=find(nPat>0);
t.extNorm2=toc; time2str(toc);
clear movNorm;


% Make Descriptor
%--------------------
nFrame=length(normPat);
clear desc
tic;
desc=cell(1,nFrame);
feaCount=zeros(1,nFrame);
for ii=1:nFrame
  if size(normPat{ii},3)>1
    feaCount(ii)=size(normPat{ii},3);
    [desc{ii},spinwavParm]=spinwavDesc(normPat{ii});
  else
    wdesc{ii}=[];
    wcoef{ii}=[];
  end
end
t.spinwavDesc=toc; time2str(toc);
clear normPat
descMat=cell2mat(desc);
[desc1,desNorm]=norm1ize(descMat);
%[wdescMat,wStat]=standize(wdescMat);



% Cluster Descriptors
%--------------------
nclust=50;
tic;
%memb=spectclust(desc1(:,1:1000),nclust);
memb=kmeanclust(desc1,nclust);
t.clust=toc;time2str(toc);
feaDist=hist(memb,nclust);
bar(feaDist,'y');axis tight;

[dummy,maxClust]=max(feaDist);
clustNum=maxClust
feaClust=descMat(:,find(memb==clustNum));
feaIm=spinwavInv(feaClust,spinwavParm);
playmov(feaIm);






% Build histogram via ssim metric
%--------------------
tic
[whist,feav]=wdescHist(wcoef,rmLev2,0.7);
nfea=length(feav);
feaim=zeros(psiz,psiz,nfea);
for ii=1:nfea
  feaim(:,:,ii)=cwtInv(feav{ii});
end
t.hist=toc; time2str(toc);
showFrame(feaim);
%playmov(feaim,0.7);





% Curvelet Regression
%--------------------
load dat.mat; 
clear optStim;
optStim.nLev=4;
optStim.nAngle=8;
optStim.useWavelet=0;
optStim.NLfeatureMap=1;
optStim.NLdimExpand=1;
optStim.NLcompressive=2
nlag=4;

[dat,stimInvOpt]=preProcStim(dat,optStim);
dat=preProcAlign(dat,0);

eStep=std(dat2.eResp)/2000;
opt.maxStep=100000;
opt.showStep=0;
opt.backTrack=100;

tic;
[bb,nIter,opt]=lagBoosting(dat.eStim,dat.eResp,eStep,dat.cStim,dat.cResp,opt);
t.boost=toc;time2str(t.boost); 

stimInvOpt.nLag=nlag;
[r,pResp,sIm]=linearizedPred(bb,dat.vStim,dat.vResp,stimInvOpt);
r
showFrame(sIm,0,1);






% Do regression
%--------------------
nlag=10
X=reshape([whist;whist],76,8000);
%X=tapdelay(X,nlag)';
X(:,1:nlag-1)=[];
Y=resp(1:end-nlag+1);
coef=X'*Y;


[bb,nIter]=l2boosting(X(1:7000,:),Y(1:7000),0.001,X(7001:end,:),Y(7001:end));




% Reconstruct Feature movies
%--------------------
tic;
if feature==1
  [feaMov1,coefCel]=spinwavInv(desc,spinwavParm);
elseif feature==2
  [feaMov2,coefCel]=spinwavInv(desc,spinwavParm);
end
t.spinwavInv=toc; time2str(toc);
%playmovs(0.5,feaMov1,feaMov2);
clear feaMov1 feaMov2 coefCel




ncmp=10;
a=randperm(length(coefCel));
[ssimMat,im]=pairssim(coefCel(a(1:ncmp)),rmLev2);
dispPairCmp(im,ssimMat);

[angMat]=zeros(ncmp);
for ii=1:ncmp
  for jj=ii+1:ncmp
    angMat(ii,jj)=anglebtw(im{ii},im{jj});
  end  % for jj
end  % for ii
dispPairCmp(im,angMat./180));















% Test curvelet decomposition & reconstruction
tic;
movR=zeros(nDim,nDim,nFrame);
for ii=1:nFrame
 ccoef=crvLet(mov(:,:,ii),0,2,8);
 movR(:,:,ii)=crvInv(ccoef,nDim,nDim);
end
figure;plot(mean(reshape(mov-abs(movR),nDim*nDim,nFrame)));
time2str(toc);



frameIdx=92;
testFrame=mov(:,:,frameIdx);
imagesc(testFrame);colormap gray;


tic
for ii=1:16
  %ii/10+0.4
  scaleFac=ii/10;
  frameResize{ii}=imresize(imresize(testFrame,scaleFac,'bicubic'), ...
    [nDim,nDim],'bicubic');
  ccoef{ii}=crvLet(frameResize{ii},0,2,8);
end
for ii=1:11
  figure(1);
  imagesc(frameResize{ii});colormap gray;
  title(ii);

  figure(2);
  crvScale(ccoef{ii},2);
  drawnow;
  pause;

end

time2str(toc);









if 0

+++++++++++++++++++++++++++++++++++++++++++++++
img = double(imread('Lena.jpg'));
n = size(img,1);
sigma = 20;        
is_real = 1;

noisy_img = img + sigma*randn(n);

disp('Compute all thresholds');
F = ones(nDim);
X = fftshift(ifft2(F)) * sqrt(prod(size(F)));
tic, C = fdct_wrapping(X,0,2); toc;

% Compute norm of curvelets (exact)
E = cell(size(C));
for s=1:length(C)
  E{s} = cell(size(C{s}));
  for w=1:length(C{s})
    A = C{s}{w};
    E{s}{w} = sqrt(sum(sum(A.*conj(A))) / prod(size(A)));
  end
end

% Take curvelet transform
disp(' ');
disp('Take curvelet transform: fdct_wrapping');
tic; C = fdct_wrapping(noisy_img,1,2); toc;

% Apply thresholding
Ct = C;
for s = 2:length(C)
  thresh = 3*sigma + sigma*(s == length(C));
  for w = 1:length(C{s})
    Ct{s}{w} = C{s}{w}.* (abs(C{s}{w}) > thresh*E{s}{w});
  end
end


end  % if 0

