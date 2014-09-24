function opt=defaultOpt(optIn,defVal,defRange)
%function opt=defaultOpt(optIn,defVal,defRange)
%
% This takes in a option flag structure [optIn] and the default
% values of each of the fields [defVal]. It will set the option
% fields to the default value when the field in [optIn] does not
% exist, or is not specified or empty or not of the same class
% as the default value, or not in the range specified in [defRange].
% Otherwise, the field value of the [optIn] will be used.
%
% EXAMPLE:
%     defVal.maxIter=1e5;      dRange.maxIter=[1,1e8];
%     defVal.numStep=2;        dRange.numStep=[2:2:10];
%     defVal.stepSiz=0.01;     dRange.stepSiz=[1e-8,1];
%     defVal.verbose=1;        dRange.verbose=[0,1];
%     defVal.func='mean';      dRange.func={'mean','mode','median'};
%     defVal.func2='exp';     dRange.func2={'exp','log','sin','cos'};
%     clear opt;        % Unspecified fields will be set to default 
%     opt.numStep=100;  % Out of range, will be set to default value
%     opt.stepSiz=5;    % Not in the list of allow values
%     opt.verbose=0;    % In range and OK. maxIter is not specified
%     opt.func='log2'   % Not in allowed list.
%     opt=defaultOpt(opt,defVal,dRange)
% The output should be:
%     opt =
%         maxIter: 100000
%         numStep: 2
%         stepSiz: 0.0100
%         verbose: 0
%            func: 'mean'
%           func2: 'exp'
%
% SEE ALSO: structCmp, consistentField, structUnion
%
% by Michael Wu - waftingp@uclink4.berkeley.edu (Mar 2003)
%
% ====================


opt=optIn;
fname=fieldnames(defVal);
nField=length(fname);


if class(optIn)~=class(defVal)
  opt=defVal;
else
  for ii=1:nField
    default=getfield(defVal,fname{ii});
    if isfield(optIn,fname{ii});
      fieldVal=getfield(optIn,fname{ii});
      if ~isempty(fieldVal);
        if any(strmatch(class(default),class(fieldVal),'exact')) | (isnumeric(default) & isnumeric(fieldVal));

          if nargin>2
            if isfield(defRange,fname{ii})
              rangeVal=getfield(defRange,fname{ii});
              if isnumeric(fieldVal)
                if ismember(fieldVal,rangeVal) | (fieldVal<=max(rangeVal) & fieldVal>=min(rangeVal))
                  opt=setfield(opt,fname{ii},fieldVal);
                else
                  warning(sprintf(['Field [.',fname{ii}, ...
                    '] NOT valid!\n----- Set to Default: .' ...
                    ,fname{ii},' = %g\n'],default));
                  opt=setfield(opt,fname{ii},default);
                end % if ismember
              
              elseif ischar(fieldVal)
                if ismember(fieldVal,rangeVal) 
                  opt=setfield(opt,fname{ii},fieldVal);
                else
                  warning(sprintf(['Field [.',fname{ii}, ...
                    '] NOT valid!\n----- Set to Default: .' ...
                    ,fname{ii},' = %s\n'],default));
                  opt=setfield(opt,fname{ii},default);
                end % if ismember
                 
              end  % if isnumeric
            %else
            %  opt=setfield(opt,fname{ii},fieldVal);
            end % isfield
          end % if nargin
        
        else
          warning(sprintf(['Field [.',fname{ii}, ...
            '] NOT valid!\n----- Set to Default: .' ...
			,fname{ii},' = %s\n'],default));

          opt=setfield(opt,fname{ii},default);
        end % if strmatch
      else
        opt=setfield(opt,fname{ii},default);
      end % if ~isempty
    else
      opt=setfield(opt,fname{ii},default);
    end % if isfield
  end % for ii
end % if


if nField==length(fieldnames(opt))
  opt=orderfields(opt,defVal);
end

