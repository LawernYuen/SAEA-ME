function [weight, neighbour] = init_weights(popsize, niche, objDim)
    % init_weights function initialize a pupulation of subproblems structure
    % with the generated decomposition weight and the neighbourhood
    % relationship.
    
    
    weight = initweight(objDim, popsize);
    weight = weight';
    %weight = normr(weight); % Normalize the weight vectors

    %Set up the neighbourhood.
    leng           = size(weight, 1);
    distanceMatrix = zeros(leng, leng);
    neighbour      = zeros(leng, niche);
    
    for i = 1 : leng
        for j = i + 1 : leng
            A                    = weight(i, :);
            B                    = weight(j, :);
            distanceMatrix(i, j) = (A - B) * (A - B)';
            distanceMatrix(j, i) = distanceMatrix(i, j);
        end
        [~, sindex]     = sort(distanceMatrix(i, :));
        neighbour(i, :) = sindex(1 : niche)';
    end
end