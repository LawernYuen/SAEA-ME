function pop = crowdingdistance(pop, F, nObj)

    nF = numel(F);
    for i = 1 : nF
        [~, costs]=getV(pop(F{i}), size(pop(1).Position,2), nObj);
        n = numel(F{i});
        d = zeros(n, nObj);
        for j = 1 : nObj
            [cj, cidx] = sort(costs(:, j));
            d(cidx(1), j) = inf;
            for k = 2 : n-1
                d(cidx(k), j)=abs(cj(k+1) - cj(k-1)) / abs(cj(1) - cj(end));
            end
            d(cidx(end), j) = inf;
        end
        for k = 1 : n
            pop(F{i}(k)).CrowdingDistance = sum(d(k, :));
        end
    end
end