function [parent_pop, parent_objs] = moead(problem, model_hyp, db_pop, db_objs)
    % This is the main procedure of MOEA/D
    
    %global variable definition.
    global params idealpoint itrCounter nFEs;
    
    %% Set the algorithms parameters.
    %Set up the initial setting for the MOEA/D.
    objDim = 2*problem.od;
    idealp = ones(1, objDim) * inf;
    params.F        = 0.5;
    params.CR       = 0.5;
    
    %the default values for the parameters.
    params.nr        = 2;
    params.niche     = 5;    
    params.delta     = 0.9;
    params.dmethod   = 'i_te';    

    [subproblems, neighbour] = init_weights(params.popsize, params.niche, objDim);
    %initial the subproblem's initital state.
    parent_pop = randompoint(problem, params.popsize, 'random');
    parent_objs = getObj(problem, model_hyp, db_pop, db_objs, parent_pop);
    assert(isreal(parent_objs), 'FAKE NEWS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    % Find the ideal point
    idealpoint = min(idealp, min(parent_objs));

    %% Main precedure
    itrCounter = 1;
    while ~terminate(itrCounter)
        [parent_pop, parent_objs] = evolve(parent_pop, parent_objs, subproblems, neighbour, problem, model_hyp, db_pop, db_objs);

        fprintf('Dimension: %d :: FE: %d :: iteration %d finished\n', problem.pd, nFEs, itrCounter);
        itrCounter = itrCounter + 1;
    end
end

function [parent_pop, parent_objs] = evolve(parent_pop, parent_objs, subproblems, neighbour, problem, model_hyp, db_pop, db_objs)
    global idealpoint params;
    
    for i = 1 : params.popsize
        
        if rand < params.delta
            matingindex = neighbour(i, :);
        else
            matingindex = 1 : params.popsize;
        end
            
        % New point generation using genetic operations
        ind = genetic_op(parent_pop, i, problem.domain, params, matingindex);
        obj = getObj(problem, model_hyp, db_pop, db_objs, ind);
        % Update the ideal point
        idealpoint = min(idealpoint, obj);
        
        % Update neighbours
        [parent_pop, parent_objs] = update(parent_pop, parent_objs, subproblems, matingindex, ind, obj, params, idealpoint);
        
        clear ind obj matingindex;
    end
end

function [parent_pop, parent_objs] = update(parent_pop, parent_objs, subproblems, matingindex, ind, ind_obj, params, idealpoint)
    newobj   = subobjective(subproblems(matingindex, :), ind_obj, idealpoint, params.dmethod);
    old_objs = parent_objs(matingindex, :);
    oldobj   = subobjective(subproblems(matingindex, :), old_objs, idealpoint, params.dmethod);
    
    C = newobj < oldobj;
    counter = sum(C);
    newC    = false(size(C));

    if counter <= params.nr
        newC = C;
                
        parent_pop(matingindex(:, newC), :)  = ind(ones(counter, 1), :);
        parent_objs(matingindex(:, newC), :) = ind_obj(ones(counter, 1), :);
    else
        nonzero_ind              = find(C);
        temp                     = randperm(counter); 
        nrInd                    = temp(1 : params.nr);
        newC(nonzero_ind(nrInd)) = 1;
                
        parent_pop(matingindex(:, newC), :)  = ind(ones(params.nr, 1), :);
        parent_objs(matingindex(:, newC), :) = ind_obj(ones(params.nr, 1), :);
    end
    
    clear C newC temp nrInd
end

function y = terminate(itrcounter)
    global params;
    y = itrcounter > params.maxItr;
end

function obj = CostFunction(problem, pop, model_hyp, db_pop, db_objs)

    global params meanfunc covfunc likfunc infe;

    switch params.model
        case 'vanilla'
            obj = problem.func(pop);
        case 'gp'
            [mean, s] = gp(model_hyp, infe, meanfunc, covfunc, likfunc, db_pop, db_objs, pop);
            obj  = [mean, mean-sqrt(s)];
        case 'kriging'
            if size(pop, 1) > 1
                [mean, s] = predictor(pop, model_hyp);
                obj  = [mean, s];
            else
                [mean, ~, s, ~] = predictor(pop, model_hyp);
                obj  = [mean, s];
            end
        case 'addGP'
            add_func = model_hyp{1};
            [mean, s] = add_func(pop);
            obj = [mean, s];
        case 'rbfn'
            obj = model_hyp(pop')';
        otherwise
            error('Undefined surrogate methods!')
    end    
end

function obj = getObj(problem, model, db_pop, db_objs, pop)
    global params;
    obji = cell(1, problem.od);
    obj = [];
    for i = 1:problem.od
        obji{i} = CostFunction(problem, pop(:,params.seps{i}), model{i}, ...
            db_pop(:,params.seps{i}), db_objs(:,i));
        obj = [obj, obji{i}];
    end

end

