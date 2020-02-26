function [seps, pop, obj] = group(problem)

    global nFEs;
    ind = rand(1,problem.pd).*(problem.domain(2,:)-problem.domain(1,:))+problem.domain(1,:);
    pop = repmat(ind,problem.pd,1);
    for i = 1:problem.pd
        pop(i,i) = 0;
    end
    indobj = problem.func(ind);
    nFEs = nFEs + 1;
    popobj = problem.func(pop);
    nFEs = nFEs + problem.pd;
    seps = cell(1,problem.od);
    for i = 1:problem.pd
        for j = 1:problem.od
            if indobj(j)~=popobj(i,j)
                seps{j} = [seps{j},i];
            end
        end
    end
    pop = [pop; ind];
    obj = [popobj; indobj];

end