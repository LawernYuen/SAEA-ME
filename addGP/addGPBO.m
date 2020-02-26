function [combFuncH, funcHs] = addGPBO(problem, params, decomp, db_pop, db_objs, numItr)

% Define these to avoid typos
DECOMP_KNOWN = 'known';
DECOMP_LEARN = 'learn';
DECOMP_RAND = 'random';
DECOMP_PLEARN = 'partialLearn';
DECOMP_STOCH1 = 'stoch1';

% Model building
dummyPts = zeros(0, size(db_pop,2)); % to build the GPs
MAX_THRESHOLD_EXCEEDS = 5;
NUM_ITERS_PER_PARAM_RELEARN = 25;
if isfield(decomp, 'M'), numGroups = decomp.M;
else numGroups = numel(decomp);
end
params = boProcessParams(db_pop, db_objs, params);

% The Decomposition
gpHyperParams.decompStrategy = params.decompStrategy;
if strcmp(params.decompStrategy, DECOMP_KNOWN)
    decomposition = decomp;
    numGroups = numel(decomposition);
    % Do some diagnostics on the decomposition and print them out
    relevantCoords = cell2mat(decomposition');
    numRelevantCoords = numel(relevantCoords);
    if numRelevantCoords ~= numel(unique(relevantCoords))
        error('The same coordinate cannot appear in different groups');
    end
    fprintf('# Groups: %d, %d/%d coordinates are relevant\n', ...
        numGroups, numRelevantCoords, problem.pd);
    
elseif isfield(decomp, 'M')
    % Now decomposition should have two fields d and M
    numGroups = decomp.M;
    
else % in this case the decomposition is given.
    numGroups = numel(decomp);
end

% Initialisation points
initPts = params.initPts;
initVals = params.initVals;
numInitPts = size(initPts, 1);

% Hyper parameters for the GP
gpHyperParams.decompStrategy = params.decompStrategy;
gpHyperParams.commonMeanFunc = params.commonMeanFunc;
gpHyperParams.meanFuncs = params.meanFuncs;
gpHyperParams.commonNoise = params.commonNoise;
gpHyperParams.noises = params.noises;
gpHyperParams.fixPr = params.fixPr;
gpHyperParams.useSamePr = params.useSamePr;
gpHyperParams.sigmaPrRange = params.sigmaPrRange;
gpHyperParams.sigmaPrRanges = params.sigmaPrRanges;
gpHyperParams.useSameSm = params.useSameSm;

% The Bandwidth
if params.useFixedBandWidth
    gpHyperParams.fixSm = true;
    if params.useSameSm
        gpHyperParams.sigmaSmRange = params.sigmaSmRange;
    else
        gpHyperParams.sigmaSmRanges = params.sigmaSmRanges;
    end
else % This BO algorithm will set the bw via its own procedure
    alBWLB = params.alBWLB;
    alBWUB = params.alBWUB;
    % Set an initial bw. This will change as the algorithm progresses
    alCurrBW = alBWUB;
end

% Define the following before proceeding
boQueries = [initPts; zeros(numItr, size(initPts,2))];
boVals = [initVals; zeros(numItr, 1)];
history = [max(initVals(cumsum(triu(ones(length(initVals))))))'; zeros(numItr, 1)];
threshExceededCounter = 0;
if isempty(initVals)
    currMaxVal = -inf;
    currMaxPt = [];
else
    [currMaxVal, maxIdx] = max(initVals);
    currMaxPt = initPts(maxIdx, :);
end
if strcmp(params.decompStrategy, DECOMP_STOCH1)
    groupSampleCounts = zeros(decomp.dMax, 1);
end

fprintf('Peforming BO (dim = %d)\n', problem.pd);
for boIter = 1 : numItr
    
    % Prelims
    numBoPts = numInitPts + boIter - 1;
    currBoQueries = boQueries(1:numBoPts, :);
    currBoVals = boVals(1:numBoPts);
    
    % First redefine ranges for the GP bandwidth if needed
    if ~params.useFixedBandWidth & ( mod(boIter-1, NUM_ITERS_PER_PARAM_RELEARN) == 0 | ...
            threshExceededCounter == MAX_THRESHOLD_EXCEEDS )
        
        if threshExceededCounter == MAX_THRESHOLD_EXCEEDS
            alBWUB = max(alBWLB, 0.9 * alCurrBW);
            threshExceededCounter = 0;
            fprintf('Threshold Exceeded %d times - Reducing BW\n', MAX_THRESHOLD_EXCEEDS);
        else
            alBWUB = max(alBWLB, 0.9 * alBWUB);
        end
        
        % Define the BW range for addGPMargLikelihood
        if alBWUB == alBWLB
            gpHyperParams.fixSm = true;
            gpHyperParams.sigmaSmRanges = alBWLB * ones(numGroups, 1);
        else
            gpHyperParams.fixSm = false;
            % Use same bandwidth for now. TODO: modify to allow different bandwidths
            gpHyperParams.useSameSm = true;
            gpHyperParams.sigmaSmRange = [alBWLB, alBWUB];
        end % alBWUB == alBWLB
        
        % Obtain the optimal GP parameters
        if ~strcmp(params.decompStrategy, DECOMP_STOCH1)
            [~, ~, ~, ~, ~, ~, ~, alCurrBWs, alCurrScales, ~, learnedDecomp, margLikeVal] = ...
                addGPDecompMargLikelihood( currBoQueries, currBoVals, dummyPts, ...
                decomp, gpHyperParams );
            
            alCurrBW = alCurrBWs(1); %TODO: modify to allow different bandwidths
            if problem.pd < 24
                learnedDecompStr = '';
                for i = 1:numel(learnedDecomp)
                    learnedDecompStr = [learnedDecompStr mat2str(learnedDecomp{i})];
                end
            else
%                 learnedDecompStr = sprintf('d<=%d', max(decomp));
                learnedDecompStr = '';
                for i = 1:numel(learnedDecomp)
                    learnedDecompStr = [learnedDecompStr mat2str(learnedDecomp{i})];
                end
            end
            fprintf(['Picked bw: %0.4f (%0.2f, %0.2f), Scale: %0.4f, ML: %0.4f, ', ...
                'Decomp: %s (%d)\n'], ...
                alCurrBW, alBWLB, alBWUB, alCurrScales(1), margLikeVal/numBoPts, ...
                learnedDecompStr, numel(learnedDecomp));
            
        else
            groupBWs = cell(decomp.dMax, 1);
            groupScales = cell(decomp.dMax, 1);
            groupDecomps = cell(decomp.dMax, 1);
            groupMargLikeVals = zeros(decomp.dMax, 1);
            gpHyperParams.decompStrategy = DECOMP_PLEARN;
            
            for d = 1:decomp.dMax
                numFullGroups = floor(problem.pd/d);
                currDecomp = d * ones(numFullGroups, 1);
                numRemDims = problem.pd - numFullGroups * d;
                if numRemDims > 0
                    currDecomp = [currDecomp; numRemDims];
                end
                
                [~, ~, ~, ~, ~, ~, ~, bws, scales, ~, currDecomp, margLikeVal] = ...
                    addGPDecompMargLikelihood(currBoQueries, currBoVals, dummyPts, ...
                    currDecomp, gpHyperParams);
                groupBWs{d} = bws;
                groupScales{d} = scales;
                groupDecomps{d} = currDecomp;
                groupMargLikeVals(d) = margLikeVal; % / size(currBoVals, 1);
                
            end
            groupSampleProbs = exp(groupMargLikeVals - max(groupMargLikeVals));
            groupSampleProbs = 0.80 * groupSampleProbs/sum(groupSampleProbs) + ...
                0.20 * ones(decomp.dMax,1)/decomp.dMax;
            fprintf('Sampling Probs: %s, Counts: %s\n', ...
                mat2str(round(groupSampleProbs', 2)), ...
                mat2str(groupSampleCounts') ...
                );
            
        end % end if strcmp(params.decompStrategy, DECOMP_STOCH1)
        
        
    end % if ~params.useFixedBandWidth ...
    
    % If stochastic pick a current GP
    if ~strcmp(params.decompStrategy, DECOMP_STOCH1)
        currGPParams.sigmaSms = alCurrBWs;
        currGPParams.sigmaPrs = alCurrScales;
        currIterDecomp = learnedDecomp;
    else
        % pick a random index.
        dIdx = sampleFromMultinomial(groupSampleProbs);
        groupSampleCounts(dIdx) = groupSampleCounts(dIdx) + 1;
        currGPParams.sigmaSms = groupBWs{dIdx};
        currGPParams.sigmaPrs = groupScales{dIdx};
        currIterDecomp = groupDecomps{dIdx};
    end
    
    % Now build the GP
    currGPParams.commonMeanFunc = gpHyperParams.commonMeanFunc;
    currGPParams.meanFuncs = gpHyperParams.meanFuncs;
    currGPParams.commonNoise = gpHyperParams.commonNoise;
    currGPParams.noises = gpHyperParams.noises;
    [~,~,~,~, combFuncH, funcHs] = addGPRegression(currBoQueries, ...
        currBoVals, dummyPts, currIterDecomp, currGPParams);
        
end % end for boIter

end

