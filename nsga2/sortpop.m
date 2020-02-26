function [pop, F] = sortpop(pop)

    % Sort Based on Crowding Distance
    [~, cdidx] = sort([pop.CrowdingDistance], 'descend');
    pop = pop(cdidx);
    
    % Sort Based on Rank
    [~, ridx] = sort([pop.Rank]);
    pop = pop(ridx);
    
    % Update Fronts
    rank = [pop.Rank];
    maxrank = max(rank);
    F = cell(maxrank, 1);
    for i = 1 : maxrank
        F{i} = find(rank == i);
    end

end