function [best_ind, best_ind_obj] = de(problem, model_hyp, db_pop, db_objs, f_min)

    global nFEs params meanfunc covfunc likfunc inf;

    %% initialisation
    pop = randompoint(problem, params.popsize, 'random');
    switch params.model
        case 'vanilla'
            pop_obj = problem.func(pop);
        case 'gp'
            [obj, s] = gp(model_hyp, inf, meanfunc, covfunc, likfunc, db_pop, db_objs, pop);
            pop_obj  = acquisition_func(f_min, obj, s);
        case 'kriging'
            [obj, s] = predictor(pop, model_hyp);
            pop_obj  = acquisition_func(f_min, obj, s);
        otherwise
            error('Undefined surrogate methods!')
    end
    
    % find the solution having the best acquisition function
    [best_ind_obj, best_idx] = min(pop_obj);
    best_ind = pop(best_idx, :);
    
    %% DE procedure (DE/rand/1)
    for itr = 2 : params.maxItr
%         comp_ind = zeros(params.popsize, 5);
%         comp_s = zeros(params.popsize, 1);
%         comp_obj = zeros(params.popsize, 1);
%         comp_pobj = zeros(params.popsize, 1);
        for i = 1 : params.popsize
            oldpoint = pop(i, :);
            idx = randperm(params.popsize);
            idx(idx == i) = [];
            a = idx(1);
            b = idx(2);
            c = idx(3);
           
            % mutation
            newpoint = pop(a, :) + params.F .* (pop(b, :) - pop(c, :));
            % crossover
            jrandom = ceil(rand * problem.pd);
            randomarray         = rand(problem.pd, 1);
            deselect            = randomarray < params.CR;
            deselect(jrandom)   = true;
            newpoint(~deselect) = oldpoint(~deselect);
           
            % repair the solution
            newpoint = max(newpoint, problem.domain(1, :));
            newpoint = min(newpoint, problem.domain(2, :));
           
            switch params.model
                case 'vanilla'
                    newpoint_obj = problem.func(newpoint);  % use true obj evaluation
                case 'gp'
                    [newobj, news] = gp(model_hyp, inf, meanfunc, covfunc, likfunc, db_pop, db_objs, newpoint);
                    newpoint_obj   = acquisition_func(f_min, newobj, news);
%                     comp_ind(i,:) = newpoint;
%                     comp_s(i,:) = news;
%                     comp_obj(i,:) = problem.func(newpoint);
%                     comp_pobj(i,:) = newobj;
                case 'kriging'
                    [newobj, ~, news, ~] = predictor(newpoint, model_hyp);
                    newpoint_obj = acquisition_func(f_min, newobj, news);
                otherwise
                    error('Undefined model type!')
            end
                
            if newpoint_obj < pop_obj(i, :)
                pop(i, :)     = newpoint;
                pop_obj(i, :) = newpoint_obj;
                
                if pop_obj(i, :) < best_ind_obj
                    best_ind     = pop(i, :);
                    best_ind_obj = pop_obj(i, :);
                end
            end
        end
        fprintf('FEs = %d :: DE_itr = %d :: Best EI = %f :: Best Fitness = %f\n', nFEs, itr, best_ind_obj, f_min);
%         clf;
%         area(pop_obj);
%         pause(0.01);
    end
    clear idx pop pop_obj s obj newpoint oldpoint randomarray db_pop db_objs
    pack
end