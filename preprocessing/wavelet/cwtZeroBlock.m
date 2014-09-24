function [cfs]=cwtZeroBlock(cfs,lev,ori,c)
%function [cfs]=cwtZeroBlock(cfs,lev,ori,c)
%
% Zero or set to [c] several blocks of CWT coeffs in [cfs] 
% specified by the scales in [lev] and orientations in [ori].
% There must be a 1-to-1 correspondance between [lev] and 
% [ori]. If not output arguments are specified, it will 
% display the coeffs subimages.
%
% INPUT:
%  [cfs] = cell array of coefficients at each scale level, 
%          with orientated coeff of each level is stored as 
%          a structure array inside the cell array of scales.
%  [lev] = Array of levels to zero out.
%  [ori] = Array of orientations at scale [lev] to zero out. 
%          It can be a array of angles in degrees. Valid 
%          orientation are +/-15, +/-45 and +/-75 degrees. If 
%          [ori] is not given at these value, it is round to 
%          the nearest valid orientation. [ori] can also be a
%          char array of  
%             'v' = (+/-)75 deg, 
%             'd' = (+/-)45 deg, 
%             'h' = (+/-)75 deg 
%          or 'l' = lowpass.
% OUTPUT:
%  [cfs] = same as input [cfs], but with zero-out coeffs.
%
% SEE ALSO: cwtZeroLev, cwtZeroOri, cwtShowStr
%
% By Michael Wu  --  waftingpetal@yahoo.com (Mar 2007)
%
% ====================


% Check Input [lev]
%--------------------
if nargin<4
  c=0;
end

nzLev=length(lev);
nzOri=length(ori);
if nzLev~=nzOri
  error('zroLevOri >< [lev] and [ori] must have same length');
end

cfLev=length(cfs)-1;
totLev=cfLev+1;
if any(lev>totLev)
  error('zroLevOri >< [lev] exceed total number of levels');
end



% Check Input [ori]
%--------------------
if ischar(ori)  % If [ori] is char array
  ori=lower(ori);
  for ii=1:nzOri
    if all(lower(ori(ii))~='vhdl')
      error('zroLevOri >< [ori] must be a char of ''v'', ''h'', ''d'' & ''l'' ');
    end  % if all
  end  % for ii
  
  oriZ=[];
  lev=lev(:);
  levTemp=[lev,zeros(nzLev,1)];
  
  for ii=1:nzOri
    switch ori(ii)
      case 'h'
        oriZ=[oriZ,15,-15];
        levTemp(ii,2)=levTemp(ii,1);
      case 'd'
        oriZ=[oriZ,45,-45];
        levTemp(ii,2)=levTemp(ii,1);
      case 'v'
        oriZ=[oriZ,75,-75];
        levTemp(ii,2)=levTemp(ii,1);
      case 'l'
        if levTemp(ii,1)~=totLev
          warning('zroLevOri >< [ori]=''l'' only for highest level');
        end
        oriZ=[oriZ,nan];
        levTemp(ii,1)=totLev;
    end  % switch ori
  end  % for ii
  
  levTemp=reshape(levTemp',1,2*nzLev);
  keepIdx=find(levTemp>0);
  lev=levTemp(keepIdx);
  nzLev=length(lev);

elseif isreal(ori)  % If [ori] is numeric array
  validAng=[15,45,75,-15,-45,-75];
  cplxAng=exp(i*validAng*pi/180);
  validVec=[real(cplxAng);imag(cplxAng)];

  [ori,oriVec]=dir2ori(ori);
  oriZ=zeros(1,nzOri);
  for ii=1:nzOri
    [dummy,oriIdx]=nearVec(oriVec(:,ii),validVec);
    oriZ(ii)=validAng(oriIdx);
  end  % for ii
end  % if ischar



% Zero out a block of coefficients at [lev] and [ori]
%--------------------
for ii=1:nzLev
  if lev(ii)<totLev  % Oriented Coeffs
    switch oriZ(ii)
      case 15
        zblk=repmat(c,size(cfs{lev(ii)}.p15));
        cfs{lev(ii)}.p15=zblk;
      case 45
        zblk=repmat(c,size(cfs{lev(ii)}.p45));
        cfs{lev(ii)}.p45=zblk;
      case 75
        zblk=repmat(c,size(cfs{lev(ii)}.p75));
        cfs{lev(ii)}.p75=zblk;
      case -15
        zblk=repmat(c,size(cfs{lev(ii)}.n15));
        cfs{lev(ii)}.n15=zblk;
      case -45
        zblk=repmat(c,size(cfs{lev(ii)}.n45));
        cfs{lev(ii)}.n45=zblk;
      case -75
        zblk=repmat(c,size(cfs{lev(ii)}.n75));
        cfs{lev(ii)}.n75=zblk;
    end  % switch oriZ
  else  % Lo-Pass Coeffs
    zblk=repmat(c,size(cfs{lev(ii)}));
    cfs{lev(ii)}=zblk;
  end  % if lev
end  % for ii



% Show modified CWT coeff as subimages.
%--------------------
if nargout<1
  cwtShowStr(cfs);
end

