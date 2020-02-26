function soi = isoi(refDirects, obj, niche)

    %% Initialization
    global params;

%     niche      = 15;
    objDim     = 2;
%     refDirects = 200;
    
    % Generate reference vectors
    [weights, neighbour] = init_weights(refDirects, niche, objDim);
    
    
    % Scale
    objScaled = rescale(obj, 'InputMin', min(obj), 'InputMax', max(obj));
    
    %% Compute EMUr
    A = 1 : params.popsize;
    itr = 1;
    fronts{1} = [];
    emUs{1} = [];
    objSca = objScaled;
    while true
        emU = computeEMU(objSca, refDirects, weights);
        frontall = [fronts{:}];
        emUidx = setxor(frontall, A);
        if size(emUidx, 1) > 1
            emUidx = emUidx';
        end
        index = [emUidx, frontall];
        emUall = [emU; zeros(numel(frontall),1)];
        EMU = [emUall, index'];
        EMU = sortrows(EMU, 'descend');
        
        fronti = find(EMU(:,1));
        emUi = EMU(fronti,1);
        fronti = EMU(fronti,2);
        emUs{itr} = emUi';
        fronts{itr} = fronti';
        itr = itr + 1;
        
        frontall = [fronts{:}];
        emUidx = setxor(frontall, A);
        if size(emUidx,2)==1
            emUs{itr} = 0;
            fronts{itr} = emUidx;
            break;
        elseif isempty(emUidx)
            break;
        end
        objSca = objScaled(emUidx,:);
    end
    
%     for i = 0 : numel(fronts)-2
%         j = numel(fronts) - i;
%         emUs{j} = emUs{j} + emUs{j-1}(1);
%     end
%     for i = 1 : numel(fronts)-1
%         emUs{i} = emUs{i} + emUs{i+1}(1);
%     end
    emUr = ones(params.popsize, 1);
    emUsall = [emUs{:}];
    frontsall = [fronts{:}];
    for i = 1 : params.popsize
        emUr(i,1) = emUsall(frontsall==i);
    end
    
    %% Associate & identify & select
    assoc = ones(params.popsize, 1);
    for i = 1 : params.popsize
        frac1 = abs(weights(:,2) .* objScaled(i,1) - weights(:,1) .* objScaled(i,2));
        frac2 = sqrt(weights(:,1).^2 + weights(:,2).^2);
        distance = frac1 ./ frac2;
        [~, assoc(i,:)] = min(distance);
    end
    
    refAssoc{1} = [];
    for i = 1 : refDirects
        refAssoc{i} = find(assoc==i)';
    end
    emptyref = cellfun('isempty',refAssoc);
    delta_ne = double(emptyref');
    for i = 2 : refDirects
        delta_ne(i) = delta_ne(i-1) + delta_ne(i);
    end
    neighbour = neighbour - delta_ne;
    neighbour = neighbour(~emptyref, :);
%     assoc = (1:200)' - delta_ne;
%     assoc = flip(assoc(~emptyref));
    for i = 1 : params.popsize
        assoc(i) = assoc(i) - delta_ne(assoc(i));
    end
    refAssoc = refAssoc(~emptyref);
    neighbour = max(neighbour, 1);
    neighbour = min(neighbour, numel(refAssoc));
    S = [];
    for i = 1 : numel(refAssoc)
        [~, idx] = max(emUr(refAssoc{i},:));
        idx = refAssoc{i}(idx);
        S = [S; idx];
    end
    
    soi = [];
    for i = 1 : numel(refAssoc)
        neibor = refAssoc(neighbour(assoc(S(i)),:));
        neiidx = [neibor{:}];
        comp = emUr(S(i)) >= emUr(neiidx);
        if sum(comp) == niche
            soi = [soi; S(i)];
        end
    end

end

function emU = computeEMU(obj, refDirects, weights)

    linearU = ones(size(obj, 1), refDirects);
    emU = zeros(size(obj, 1), refDirects);
    for i = 1 : refDirects
        linearU(:,i) = weights(i,1) * obj(:,1) + weights(i,2) * obj(:,2);
        sortU = unique(linearU(:,i));
        minU = sortU(1);
        min2U = sortU(2);
        minUidx = linearU(:,i)==minU;
        emU(minUidx, i) = min2U - minU;
    end
    emU = sum(emU,2)/refDirects;

end