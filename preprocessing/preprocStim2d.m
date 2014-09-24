function [stim,params]=preprocStim2d(rawStim,params)
%function [dat,params]=preProcStim2d(dat,opt)
%
% Perform preprocessing transform to all the stimuli in [dat].
% 
% INPUT:
%           [rawStim] = A X-by-Y-by-T matrix containing stimuli (movie)
%            [params] = structure that contains parameters for preprocessing
%       .NLfeatureMap = specifies which transform to use on the raw stimuli 
%				  		1: Complex wavelet transform (Default)
%                 		2: Contourlet transform.
%                 		3: Digital curvelet transform.
%                 		4: SPin Feature transform.
% 			    .nLev = # of scale levels of wavelet, contourlet, and curvelet transform
% 	   	   .useWavlet = specifies wheter to use wavelets at the finest scale in the curvelet transform
%			 			1: USE wavelet at the finest scale = DEFAULT.
%      		         	0: DON'T use wavelet at the finest scale, use curvelet. 
%    		  .nAngle = specifies # of angles at the 2nd coarse scale in the curvelet transform 
%						[DEFAULT = minimum = 8]. Must be multiples of 4.
%  				 .pyr = scalar in [1,9] specifying the pyramidal filter type for contourlet transform
%              			1: '9-7'      
%						2: '5-3'      
%						3: 'Burt'
%              			4: 'pkva'     
%						5: 'haar'
%				     	6: 'db2'
%              			7: 'coif1'    
%						8: 'sym2'     
%						9: 'dmey' (Default)
%  				 .dir = scalar in [1,7] specifying the directional filter type for contourlet transform
%              			1: 'haar'     
%						2: '5-3' (Default)     
%              			3: '9-7'      
%						4: 'pkva'     
%						5: 'pkva6'
%              			6: 'pkva8'    
%						7: 'pkva12'
%          .phasemode = A parameter to specify how to deal with phase information
%                       0: spectral amplitude (default)
%                       1: linear sin and cos phase amplitude (2x number of channels)
%                       2: half-rectified sin and cos phase amplitude (4x number of channels)
%                       3: 0+1 (3x number of channels)
%                       4: 0+2 (5x number of channels)
% 	   	   .NLoutput = specifies output nonlinearity applied to all channels after preprocessing stimuli
%					   >0: The exponent to which all channels will be raised
%                 		0: Log Transform log(1+x)
%                	   -1: Arc Tangent
%                 		[otherwise] = No Compressive NL (Default)
%    	 .threshCoef = specifies whether to keep only the largest coefficients
%						0: Don't threshold small coeffs (Default) 
%					   >0: Other values N, will only keep N largest coeff. 
%         .normalize = specifies whethter to z-score the components of the feature vector
%						1: Z-score channels (Default)
%                 		0: Don't z-score.
%    	  .addDCterm = specifies whether to add a DC channel to the preprocessed stimuli 
%						0: Don't add a DC column of 1's (Default)
%                 		1: Add a DC column of 1's as the constant offset feature
% OUTPUT:
%             [stim] = Preprocessed stimuli that can be used for STRF fitting.
%                 		NxD matrix, N=sample size, D=dimensions.
% 			[params] = structure that contains parameters for preprocessing, with additional 
%						fields related to the chosen NL feature map and:
%              .stds = standard deviations for each channel (set if .normalize is non-zero)
%             .means = means for each channel (set if normalize is non-zero)
%    .zeromean_value = Offset value of total movie. This is set if .zeromean is non-zero.
%						
%
% SEE ALSO: preProcData.m, setStimParm.m, setRespParm.m
% 
% Modified for strfLab inclusion by Michael Oliver  -- michael.d.oliver@gmail.com (Sept 2008)
% Original function by Michael Wu  --  waftingpetal@yahoo.com (Jun 2007)
% ====================


strNLfeatureMap={ ...
  'Complex Wavelet Transform', ...
  'Contourlet Transform', ...
  'Digital Curvelet Transform' ...
  'SPin WOrD Feature Histogram'};
nFeatureMap=length(strNLfeatureMap);

strphasemode={ ...
  'Spectral Amplitude', ...
  'Linear Real and Imaginary Phase', ...
  'Half-Rectified Real and Imaginary Phase Separated', ...
  'Spectral Amplitude + Linear Real and Imaginary Phase', ...
  'Spectral Amplitude + Half-Rectified Real and Imaginary Phase Separated'};

strNLoutput={ ...
  'Power Transform ', ...
  'Log Transform log(1+x)', ...
  'Arc Tangent', ...
  'No Compressive NL'};
nCompressive=length(strNLoutput);


% Set default option values
%--------------------
% Nonlinear transform of features
optDef.NLfeatureMap=1;
optRng.NLfeatureMap=[1:nFeatureMap];
optDef.phasemode=0;
optRng.phasemode=[0:4];
optDef.NLoutput=[];
optRng.NLoutput=[-1:10];
optDef.threshCoef=0;
optRng.threshCoef=[0,inf];
optDef.normalize=1;
optRng.normalize=[0,1];
optDef.zeromean=1;
optRng.zeromean=[0,1];
optDef.addDCterm=0;
optRng.addDCterm=[0,1];
optDef.verbose=1;
optRng.verbose=[0,1];
optDef.class = 'preprocStim2d';

if nargin<1
  stim=optDef;
  return;
end
if nargin<2
  params=optDef;
else
  params=defaultOpt(params,optDef,optRng);
end

% Get Stim Size & Inits
%--------------------
if nargin>0
  movSiz=size(rawStim);
  params.frameSiz=movSiz(1:2);
end

switch params.NLfeatureMap
  case 1  % Complex Wavelet Options
    optDef.nLev=ceil(log2(min(params.frameSiz)));
    optRng.nLev=[1,16];    
  case 2
	optDef.nLev=min(find(factor(min(params.frameSiz))>2))-1;
	optRng.nLev=[1,(min(find(factor(min(params.frameSiz))>2))-1)];
	optDef.pyr=9;
	optRng.pyr=[1,9];
	optDef.dir=2;
	optRng.dir=[1,7];
  case 3  % Curvelet Transform Options
    optDef.nLev=ceil(log2(min(params.frameSiz))-3);
    optRng.nLev=[1,16];
    optDef.nAngle=8;
    optRng.nAngle=[8:4:32];
    optDef.useWavelet=0;
    optRng.useWavelet=[0,1];
end
params=defaultOpt(params,optDef,optRng);


% Feature Map
%--------------------
switch params.NLfeatureMap
  case 1
    disp(['***** ',strNLfeatureMap{params.NLfeatureMap},' *****']);
    xfunc=@mov2cwt;
    params.mat2strFunc=@cwtMat2str;
    params.stimInvFunc=@cwtInv;
  case 2
    disp(['***** ',strNLfeatureMap{params.NLfeatureMap},' *****']);
    xfunc=@mov2cont;
    params.mat2strFunc=@vec2pdfb;
    params.stimInvFunc=@contInv;
  case 3
    disp(['***** ',strNLfeatureMap{params.NLfeatureMap},' *****']);
    xfunc=@mov2crv;
    params.mat2strFunc=@crvMat2str;
    params.stimInvFunc=@crvInv;
  case 4
    disp(['***** ',strNLfeatureMap{params.NLfeatureMap},' *****']);
    xfunc=@mov2spin;
    %params.mat2strFunc=@spinMat2str;
    %params.stimInvFunc=@spinInv;    
  otherwise
    error('preProcStim2d >< Unknow Feature Transform!');
end  % switch

if params.NLfeatureMap<4
  rawStim=single(rawStim);
else
  rawStim=uint8(rawStim);
end

if params.zeromean
	if params.verbose, fprintf('***** zero mean stimuli *****\n'); end
	if isfield(params, 'zeromean_value')
		rawStim = rawStim - params.zeromean_value;
	else
		thismean = mean(rawStim(:));
		rawStim = rawStim - thismean;
		params.zeromean_value = thismean;
	end
end

[stim,params.strSiz]=feval(xfunc,rawStim,params);


% Perform Threshold of Small Coeff
%--------------------
if params.threshCoef>0
  disp(['***** Thresholding: Keeping ',num2str(params.threshCoef),' largest coeffs *****']);
  stim=rowThresh(stim,params.threshCoef);
end

params.baseChan = size(stim,2);
% Perform Nonlinear Dimensional Expansion
%--------------------
switch params.phasemode
  case 0
    disp(['***** ',strphasemode{params.phasemode+1},' *****']);
    stim= abs(stim);
  case 1
    disp(['***** ',strphasemode{params.phasemode+1},' *****']);
    if ~isreal(stim)
    	stim= [real(stim), imag(stim)];
	end
  case 2
    disp(['***** ',strphasemode{params.phasemode+1},' *****']);
    if ~isreal(stim)
		stim=max([real(stim),-real(stim),imag(stim),-imag(stim)],0);
    else
	  	stim=max([stim,-stim],0);
    end
  case 3
    disp(['***** ',strphasemode{params.phasemode+1},' *****']);
    if ~isreal(stim)
      stim=[abs(stim), real(stim),imag(stim)];
	else
	  stim=[abs(stim), stim];
    end
  case 4
    disp(['***** ',strphasemode{params.phasemode+1},' *****']);
    if ~isreal(stim)
		stim= [abs(stim) max([real(stim),-real(stim),imag(stim),-imag(stim)],0)];
    else
		stim= [abs(stim), max([stim,-stim],0)];  
    end
end  % switch


% Apply Point Compressive Nonlinearity
%--------------------
if params.NLoutput > 0
	disp(['***** ',strNLoutput{1}, num2str(params.NLoutput), ' *****']);
    stim=power(stim,params.NLoutput);
elseif params.NLoutput == 0
	disp(['***** ',strNLoutput{2},' *****']);
    stim=log(1+stim);
elseif params.NLoutput == -1
	disp(['***** ',strNLoutput{3},' *****']);
    stim=atan(stim);
else
    disp(['***** ',strNLoutput{end},' *****']);
end

% Add DC (constant) columns
%--------------------
if params.addDCterm
  disp(['***** Adding DC columns to stimulus *****']);
  neSamp=size(stim,1);
  stim=[ones(neSamp,1),stim];
end

% z-score the data
%--------------------
if params.normalize
	if isfield(params, 'means') % Already preprocessed. Use the means and stds
		global s; s = stim; clear stim; % a trick to modify stim without copying it
		norm_std_mean_global(params.stds, params.means);
		stim = s; clear s;
	else
		global s; s = stim; clear stim; % a trick to modify stim without copying it
		[stds, means] = norm_std_mean_global();
		stim = s; clear s;
		params.means = means;
		params.stds = stds;
	end
end

% Convert to Single & Store stimulus inversion parameters
%--------------------
params.nChan = size(stim,2);
stim=single(stim);
params.stimLen=size(stim,2);


