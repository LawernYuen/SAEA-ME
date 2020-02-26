function best_ind = pso(problem, model_hyp, db_pop, db_objs, f_min)

    global nFEs params meanfunc covfunc likfunc inf;

    %% initialization
    inertia = 1.0;
    correction_factor = 2.0;
    pop = randompoint(problem, params.popsize, 'random');
    pop_best = zeros(params.popsize, problem.pd);
    velocity = zeros(params.popsize, problem.pd);
    obj_best = 1000 * ones(params.popsize, 1);
    
    %% The main loop of PSO
    for itr = 1 : params.maxItr
        pop = pop + velocity/1.3;
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

        % compare the function values to find the best ones
        for i = 1 : params.popsize
            if pop_obj(i, :) < obj_best(i, :)
                pop_best(i, :) = pop(i, :);
                obj_best = pop_obj;
            end
        end

        [~, fbest] = min(obj_best);                           % find the best function value in total
        best_ind = pop(fbest, :);

        % update the velocity of the particles
        velocity = inertia * (rand(params.popsize, problem.pd) .* velocity) + correction_factor * (rand(params.popsize, problem.pd) .* (pop_best ...
            - pop)) + correction_factor * (rand(params.popsize, problem.pd) .* (pop_best(fbest, :) - pop));
        
        fprintf('FEs = %d :: Itr = %d :: Best EI = %f :: Best Fitness = %f\n', nFEs, itr, fbest, f_min);
    end
    
end