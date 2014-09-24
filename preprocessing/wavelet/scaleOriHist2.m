function [hisGram]=scaleOriHist(cfMat,strSiz)
%function [hisGram]=scaleOriHist(cfMat,strSiz)
%
% No help yet
%

% Check & Format Input
%--------------------
[cfDim,cfSamp]=size(cfMat);
strSiz=wrev(strSiz);
nLev=length(strSiz);
nOri=size(strSiz{end}.idx,1);


% Initialization
%--------------------
fineDim=4;
fineIdx{1}=[ 1  2  5  6];
fineIdx{2}=[ 3  4  7  8];
fineIdx{3}=[ 9 10 13 14];
fineIdx{4}=[11 12 15 16];
fineIdx{5}=[ 6  7 10 11];
fineIdx{6}=[ 5  6  9 10];
fineIdx{7}=[ 2  3  6  7];
fineIdx{8}=[ 7  8 11 12];
fineIdx{9}=[10 11 14 15];


% Make Scale Oriented histogram
%--------------------
hisGram=cell(nLev,nOri);
for ii=1:nLev*nOri
  hisGram{ii}=single(hisGram{ii});
end
sumIdx=strSiz{1}.idx(1,1):strSiz{1}.idx(1,2);
hisGram{1,1}=sum(cfMat(sumIdx,:),1)/(length(sumIdx).^2);

for ii=2:nLev-1
  for jj=1:nOri
    sumIdx=strSiz{ii}.idx(jj,1):strSiz{ii}.idx(jj,2);
    idxLen=length(sumIdx);
    
    if idxLen>4
      for kk=1:fineDim
        sumSubIdx=sumIdx(fineIdx{kk});
        hisGram{ii,jj}(kk,:)=sum(cfMat(sumSubIdx,:),1);%.*hisGram{ii,jj}(kk,:);
      end
    elseif idxLen>1 & idxLen<=4
      hisGram{ii,jj}=sum(cfMat(sumIdx,:),1);
      hisGram{ii+1,jj}=cfMat(sumIdx,:);
    else
      hisGram{ii,jj}=cfMat(sumIdx,:);
    end
  end
end


% Convert to Matrix
%--------------------
hisGram=hisGram(:);
hisGram=cell2mat(hisGram);


