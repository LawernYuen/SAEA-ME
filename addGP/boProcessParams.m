function paramsGP = boProcessParams(db_pop, db_objs, paramsGP)
% Check for parameters and if not use default values.
    
    % Initialisation
    fprintf('Obtaining %d Points for Initialisation.\n', numel(db_objs));
    paramsGP.initPts = db_pop;
    paramsGP.initVals = db_objs;
    stdInitVals = std(paramsGP.initVals);

    % GP Hyper parameters for Regression
    % The Common Mean Function
    if ~isfield(paramsGP, 'commonMeanFunc')
        if isfield(paramsGP, 'meanFuncValue')
            paramsGP.commonMeanFunc = @(arg) paramsGP.meanFuncValue * ones(size(arg,1), 1);
        else
            paramsGP.commonMeanFunc = []; % By default will use all zeros.
        end
    end
    % The Mean function for the individual GPs
    if ~isfield(paramsGP, 'meanFuncs') || isempty(paramsGP.meanFuncs)
        paramsGP.meanFuncs = @(arg) zeros( size(arg,1), 1);
    end
    % Common noise parameters
    if ~isfield(paramsGP, 'commonNoise') || isempty(paramsGP.commonNoise)
        paramsGP.commonNoise = 0.01 * stdInitVals;
    end
    % Individual noise
    if ~isfield(paramsGP, 'noises') || isempty(paramsGP.noises)
        paramsGP.noises = 0;
    end
    
    % Scale parameters
    if ~isfield(paramsGP, 'fixPr')
        paramsGP.fixPr = false;
    end
    if ~isfield(paramsGP, 'useSamePr')
        paramsGP.useSamePr = true;
    end
    if ~isfield(paramsGP, 'sigmaPrRange')
        paramsGP.sigmaPrRange = [0.03 30] * stdInitVals;
    end
    if ~isfield(paramsGP, 'sigmaPrRanges')
        paramsGP.sigmaPrRanges = []; % use same prior by default
    end
    
    % Bandwidth parameters
    if ~isfield(paramsGP, 'useFixedBandWidth')
       paramsGP.useFixedBandWidth = false;
    end
    if ~isfield(paramsGP, 'fixSm')
      paramsGP.fixSm = false;
    end
    if ~isfield(paramsGP, 'useSameSm')
      paramsGP.useSameSm = true;
    end
    if ~isfield(paramsGP, 'alBWLB')
       paramsGP.alBWLB = 1e-5;
    end
    if ~isfield(paramsGP, 'alBWUB')
       paramsGP.alBWUB = 5;
    end

end
