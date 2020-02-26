function knee = idtfy_knee(pop, popobj)
    [~, extre_idx] = min(popobj);
    A = popobj(extre_idx, :);
    b = ones(numel(extre_idx), 1);
    x = A \ b;
    
    popsize = size(popobj, 1);
    coff = repmat(x', popsize, 1);
    frac1 = abs(sum(popobj .* coff, 2) + 1);
    frac2 = norm(x);
    dist = frac1 / frac2;
    [~, idx] = max(dist);
    knee = pop(idx, :);
end