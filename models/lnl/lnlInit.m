%% Initialize a decoupled linear-nonlinear model. Two separate models need
%  to be specified - one a purely linear model, and another nonlinear model
%  which takes a scalar input and gives a scalar output.
%
% Input:
%
%   linModel: the linear model
%
%   nlModel: the nonlinear model
%
function modelParams = lnlInit(linModel, nlModel)

    modelParams = struct;
    modelParams.type = 'lnl';
    modelParams.nWts = linModel.nWts + nlModel.nWts;
    modelParams.linModel = linModel;
    modelParams.nlModel = nlModel;
    