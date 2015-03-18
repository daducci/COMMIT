%
% Set the default parameters for each model
%
% Parameters
% ----------
% model : string
%   String identifying the model to use
%
function COMMIT_SetModel( model )
    global CONFIG

    % Call the specific model constructor
    modelClass = str2func( ['COMMIT_' upper(model)] );
    if exist([func2str(modelClass) '.m'],'file')
        CONFIG.model = modelClass();
    else
        error( '[COMMIT_SetModel] Model "%s" not recognized', model )
    end
end
