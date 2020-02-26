function [pos, obj] = getV(pop, d, nObj)

    pos = ones(numel(pop), d);
    obj = ones(numel(pop), nObj);
    for i = 1 : numel(pop)
        pos(i,:) = pop(i).Position;
        obj(i,:) = pop(i).Cost;
    end
    
end
