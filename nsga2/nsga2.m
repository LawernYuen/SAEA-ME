function [popp, popo] = nsga2(problem, model_hyp, db_pop, db_objs)

    %global variable definition.
    global params nFEs;
    
    %% Set the algorithms parameters.
    nObj       = 2*problem.od;
    nVar       = problem.pd;
    lb         = problem.domain(1,:);
    ub         = problem.domain(2,:);
    nPop       = params.popsize;
    
    % SBX & polynomial mutation parameters
    proC = 1;    % crossover probability
    disC = 20;   % distribution index of sbx
    proM = 1;    % expectation of number of bits doing mutation
    disM = 20;   % distribution index of polynomial mutation
    
    % BLX & normally distributed mutation parameters
    pCrossover = 0.7;                         % Crossover Percentage
    nCrossover = 2*round(pCrossover*nPop/2);  % Number of Parnets (Offsprings)
    pMutation  = 0.4;                         % Mutation Percentage
    nMutation  = round(pMutation*nPop);       % Number of Mutants
    mu         = 0.02;                        % Mutation Rate
    sigma      = 0.1*(ub-lb);                 % Mutation Step Size
    
    % DE operators parameters
    deF  = 0.5;
    deCR = 0.5;

    %% Initialization
    pop  = getPop(nPop);
    popp = rand(nPop, nVar) .* (ub-lb);
    rng('shuffle');
    popo = getObj(problem, model_hyp, db_pop, db_objs, popp);
    assert(isreal(popo), 'FAKE NEWS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    pop = assignV(popp, popo, pop);
    
    [pop, F] = nonDominatedSort(pop);  % Non-Dominated Sorting
    pop = crowdingdistance(pop, F, nObj);    % Calculate Crowding Distance
    [pop, ~] = sortpop(pop);           % Sort Population
    
    %% NSGA-II Main Loop
    itr = 0;
    F1count = 0;
    while true
        itr = itr + 1;
        popPos = getV(pop, problem.pd, nObj);

        %% (DE/rand/1)
        idx = randperm(params.popsize);
        idx(idx == 1) = [];
        for i = 2 : params.popsize
            idxtemp = randperm(params.popsize);
            idxtemp(idxtemp == i) = [];
            idx = [idx; idxtemp];
        end
        a = idx(:, 1);
        b = idx(:, 2);
        c = idx(:, 3);

        % Mutation
        newpoint = popPos(a,:) + deF*(popPos(b,:) - popPos(c,:));
        % Crossover
        jrandom             = ceil(rand(params.popsize,1) * problem.pd);
        randomarray         = rand(params.popsize, problem.pd);
        deselect            = randomarray < deCR;
        linearInd           = sub2ind(size(deselect),1:params.popsize,jrandom');
        deselect(linearInd) = true;
        newpoint(~deselect) = popPos(~deselect);

        % repair
        newpoint = max(newpoint, lb);
        newpoint = min(newpoint, ub);

        newobj = getObj(problem, model_hyp, db_pop, db_objs, newpoint);
        popnew = getPop(params.popsize);
        popnew = assignV(newpoint, newobj, popnew);
        rng('shuffle');

        %% BLX crossover
        popc = getPop(nCrossover);
        i1 = randi(nPop, [nCrossover/2, 1]);
        i2 = randi(nPop, [nCrossover/2, 1]);
        p1 = popPos(i1, :);
        p2 = popPos(i2, :);
        [popcPos1, popcPos2] = crossover(p1, p2, lb, ub);
        popcPos = [popcPos1; popcPos2];
        popcObj = getObj(problem, model_hyp, db_pop, db_objs, popcPos);
        popc = assignV(popcPos, popcObj, popc);

        % Mutation        
        popm = getPop(nMutation);
        i = randi(nPop, [nMutation, 1]);
        p = popPos(i, :);
        popmPos = mutate(p, mu, sigma, lb, ub);
        popmObj = getObj(problem, model_hyp, db_pop, db_objs, popmPos);
        popm = assignV(popmPos, popmObj, popm);
        rng('shuffle');

        %% SBX & polynomial mutation
        offspring = reproduction(popPos, proC, disC, proM, disM, lb, ub);
        offspringCost = getObj(problem, model_hyp, db_pop, db_objs, offspring);
        popr = getPop(size(offspring, 1));
        popr = assignV(offspring, offspringCost, popr);
        rng('shuffle');

        %%

        % Merge
        pop = [pop; popnew; popc; popm; popr];
        [pop, F] = nonDominatedSort(pop);  % Non-Dominated Sorting
        pop = crowdingdistance(pop, F, nObj);    % Calculate Crowding Distance
        pop = sortpop(pop);                % Sort Population
        pop = pop(1:nPop);                 % Truncate

        [pop, F] = nonDominatedSort(pop);  % Non-Dominated Sorting
        pop = crowdingdistance(pop, F, nObj);    % Calculate Crowding Distance`
        pop = sortpop(pop);                % Sort Population
        F1 = pop(F{1});                    % Store F1
        
        disp(['Dimension = ' num2str(problem.pd) ': FE = ' num2str(nFEs)...
            ': Iteration ' num2str(itr) ': Number of F1 Members = ' num2str(numel(F1))]);
        
        % Counting stopping criteria
        if numel(F{1}) == params.popsize
            F1count = F1count + 1;
        else
            F1count = 0;
        end
        
        % Stopping criteria
        if (F1count >= 20 && itr >= params.maxItr) || itr >= 300
            break;
        end
    end

    [popp, popo] = getV(F1, problem.pd, nObj);
    
    clear a b c db_pop db_objs deCR deF deselect disC disM F F1 F1count i i1 i2 idx idxtemp itr ...
        jrandom lb linearInd model_hyp mu nCrossover newobj newpoint nMutation nObj nPop nVar ...
        offspring offspringCost p p1 p2 pCrossover pMutation pop popc popcObj popcPos popcPos1 ...
        popcPos2 popm popmObj popmPos popnew popPos popr problem proC proM randomarray ...
        ub
    
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
                obj  = [mean, mean-sqrt(s)];
            else
                [mean, ~, s, ~] = predictor(pop, model_hyp);
                obj  = [mean, mean-sqrt(s)];
            end
        case 'addGP'
            add_func = model_hyp{1};
            [mean, s] = add_func(pop);
            obj = [mean, mean-sqrt(s)];
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
