function soi = isoiV2(refDirects, obj, niche)

    %% Initialization

%     niche      = 15;
    objDim     = 2;
    nPop = size(obj, 1);
%     refDirects = 200;
    
    % Generate reference vectors
    [weights, neighbour] = init_weights(refDirects, niche, objDim);
    
    
    % Scale
    objScaled = rescale(obj, 'InputMin', min(obj), 'InputMax', max(obj));
    
    %% Compute EMUr
    itr = 1;
    emuMax= 0;
    empty_individual.obj = [];          % objective value
    empty_individual.objScaled = [];    % scaled objective value
    empty_individual.emu = [];          % emu
    empty_individual.emur = [];         % emur
    empty_individual.front = [];        % emu front
    empty_individual.ref = [];          % associated weight vector
    pop = repmat(empty_individual, nPop, 1);
    
    for i = 1:nPop
        pop(i).obj = obj(i,:);
        pop(i).objScaled = objScaled(i,:);
        pop(i).emu = 0;
    end
    while true
        popEmu = [];
        for i = 1:nPop
            popEmu = [popEmu; pop(i).emu];
        end
        index = find(~popEmu);
        if numel(index)==0
            break;
        elseif numel(index)==1
            pop(index).front = itr;
            pop(index).emur = emuMax(itr);
            break;
        end
        emu = computeEMU(objScaled(index,:), refDirects, weights);
        emuMax = [emuMax; max(emu)];
        FIndex = find(emu);
        for i = 1:numel(FIndex)
            pop(index(FIndex(i))).emu = emu(FIndex(i));
            pop(index(FIndex(i))).emur = emu(FIndex(i)) + emuMax(itr);
            pop(index(FIndex(i))).front = itr;
        end
        disp(['Iteration of computing EMUr ' num2str(itr) ' finished']);
        itr = itr + 1;
    end
    
    %% Associate & identify
    solAssoc = cell(refDirects, 1);         % solutions associated to each weight vector
    solIdtf = zeros(refDirects, 1);         % solutions selected at identify step
    emur = [];                              % emur of pop
    solRefNeighbor = cell(refDirects, 1);   % solutions of neighbor weights of each weight
    soia = [];                              % SOIa, most preferred set of solutions
    for i = 1:nPop
        distance = calDist(weights, objScaled(i,:));
        [~, refIdx] = min(distance);
        pop(i).ref = refIdx;
        solAssoc{refIdx} = [solAssoc{refIdx}, i];
        emur = [emur; pop(i).emur];
    end
    for i = 1:refDirects
        if ~isempty(solAssoc{i})
            emurIdtf = emur(solAssoc{i});
            [~, idxtemp] = max(emurIdtf);
            solIdtf(i) = solAssoc{i}(idxtemp);
        end
        for j = 2:niche
            solRefNeighbor{i} = [solRefNeighbor{i}, solAssoc{neighbour(i,j)}];
        end
    end
    
    %% Select
    nonzeroIdx = find(solIdtf);
    for i = 1:numel(nonzeroIdx)
        if emur(solIdtf(nonzeroIdx(i))) >= max(emur(solRefNeighbor{nonzeroIdx(i)}))
            soia = [soia; solIdtf(nonzeroIdx(i))];
        end
    end
    soi = soia;
    
end

function emU = computeEMU(obj, refDirects, weights)

    linearU = ones(size(obj, 1), refDirects);
    emU = zeros(size(obj, 1), refDirects);
    for i = 1 : refDirects
        linearU(:,i) = weights(i,1) * obj(:,1) + weights(i,2) * obj(:,2);
        sortU = unique(linearU(:,i));
        minU = sortU(1);
        if numel(sortU) == 1
            min2U = minU;
        else
            min2U = sortU(2);
        end
        minUidx = linearU(:,i)==minU;
        emU(minUidx, i) = min2U - minU;
    end
    emU = sum(emU,2)/refDirects;

end

function distance = calDist(weight, obj)

    frac1 = abs(weight(:,2) .* obj(1) - weight(:,1) .* obj(2));
    frac2 = sqrt(weight(:,1).^2 + weight(:,2).^2);
    distance = frac1 ./ frac2;

end