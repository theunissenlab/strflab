function [hisGram]=scaleOriHist(cfMat,strSiz)
%function [hisGram]=scaleOriHist(cfMat,strSiz)
%
% No help yet
%



% Check & Format Input
%--------------------
[cfDim,cfSamp]=size(cfMat);
strMat=cell2mat(strSiz);
hisIdx=cat(1,strMat.idx);
nOri=size(strSiz{1}.idx,1);
fineDim=4;


% Initialization
%--------------------
fineIdx{1}=[ 1  2  5  6];
fineIdx{2}=[ 3  4  7  8];
fineIdx{3}=[ 9 10 13 14];
fineIdx{4}=[11 12 15 16];
fineIdx{5}=[ 6  7 10 11];
fineIdx{6}=[ 5  6  9 10];
fineIdx{7}=[ 2  3  6  7];
fineIdx{8}=[ 7  8 11 12];
fineIdx{9}=[10 11 14 15];


% Compute Summed Index
%--------------------
sumIdx={};
for jj=1:nOri
  curIdx=hisIdx(jj,1):hisIdx(jj,2); 
  for kk=1:fineDim
    sumIdx{fineDim*(jj-1)+kk}=curIdx(fineIdx{kk});
  end
end

hisIdx(1:nOri,:)=[];
lastIdx=fineDim*nOri;

for ii=1:size(hisIdx)
  sumIdx{lastIdx+ii}=hisIdx(ii,1):hisIdx(ii,2);
end


% Make Orientation histogram
%--------------------
hisDim=length(sumIdx);
hisGram=zeros(hisDim,cfSamp);
for ii=1:hisDim
  hisGram(ii,:)=sum(cfMat(sumIdx{ii},:),1);
end
hisGram(hisDim,:)=hisGram(hisDim,:)/(length(sumIdx{hisDim}).^2);


