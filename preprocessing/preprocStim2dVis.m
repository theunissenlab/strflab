function [out, strf] = preprocStim2dVis(strf, visopt)

nDelay=length(strf.delays);
reps = strf.params.nChan / strf.params.baseChan;

strfIm.excit=zeros(strf.params.frameSiz(1),strf.params.frameSiz(2),nDelay);
strfIm.inhib=zeros(strf.params.frameSiz(1),strf.params.frameSiz(2),nDelay);

if strf.params.NLfeatureMap == 2
	varargin{2} = strf.params;
elseif strf.params.NLfeatureMap == 3
	varargin{2} = strf.params.frameSiz;
end


  for ii=1:nDelay
	for iii = 1:reps
		switch strf.params.phasemode
		  case 0
			roicvt = 1;
		  case 1
		    if iii == 2
				roicvt = 1*i;
			else
				roicvt = 1;
			end
		  case {2, 3}
		    if iii > 2
				roicvt = 1*i;
			else
				roicvt = 1;
			end
		  case 4
		    if iii > 3
				roicvt = 1*i;
			else
				roicvt = 1;
			end
		end  % switch

		cfw = strf.w1(((iii-1)*strf.params.baseChan)+1:(iii*strf.params.baseChan),1,ii);
    	cf=squeeze(cfw);
    	cfSign=cf.*(cf>0);
		cfSign = cfSign .* roicvt;
    	cfStr=feval(strf.params.mat2strFunc,cfSign,strf.params.strSiz);  % invert to coeff structure
		varargin{1} = cfStr;
		strfIm.excit(:,:,ii)= strfIm.excit(:,:,ii) + real(feval(strf.params.stimInvFunc,varargin{:}));
  
    	cfSign=cf.*(cf<0);
		cfSign = cfSign .* roicvt;
    	cfStr=feval(strf.params.mat2strFunc,cfSign,strf.params.strSiz);  % invert to coeff structure
		varargin{1} = cfStr;
    	strfIm.inhib(:,:,ii)= strfIm.inhib(:,:,ii) + real(feval(strf.params.stimInvFunc,varargin{:}));

	end
  end

strfImDisp=cat(1,strfIm.excit,strfIm.inhib);
srfSiz=size(strfImDisp);
nSrfFrm=srfSiz(3);
rshSiz=srfSiz(2)*srfSiz(3);

smax=max(strfImDisp(:));
smin=min(strfImDisp(:));
if smax==smin
  smax=smax+1;
  smin=smin-1;
end

% clf(hf);
% rfRadius=(dat.opt.stim.stimSiz/dat.parm.stimwindowcrf)/2;
% rfCenter=dat.opt.stim.stimSiz/2+0.5;

imagesc(reshape(strfImDisp,srfSiz(1),rshSiz),[-smax,smax]);
axis image; colormap gray;

set(gca, 'XTick', 0:strf.params.frameSiz(2):strf.params.frameSiz(2)*size(strf.delays,2))
set(gca, 'YTick', 0:strf.params.frameSiz(1):strf.params.frameSiz(1)*2)
set(gca, 'XTickLabel', []) 
set(gca, 'YTickLabel', [])
set(gcf, 'position', [1300 1100 1200 400])
grid on;
xlabel('time (frame)');
ylabel('inhibitory / excitatory');

