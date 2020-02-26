function best_ind = direct(problem, model_hyp, db_pop, db_objs, f_min)
    

    global nFEs params meanfunc covfunc likfunc inf;

    %% initialisation
    pop = randompoint(problem, 1, 'random');
    best_ind = patternsearch(@objfun, pop);
    fprintf('FEs = %d :: Best EI = %f :: Best Fitness = %f\n', nFEs, objfun(best_ind), f_min);
    
    function obj = objfun(x)
        switch params.model
            case 'vanilla'
                obj = problem.func(x);
            case 'gp'
                [objv, s] = gp(model_hyp, inf, meanfunc, covfunc, likfunc, db_pop, db_objs, x);
                obj  = acquisition_func(f_min, objv, s);
            case 'kriging'
                [objv, s] = predictor(x, model_hyp);
                obj  = acquisition_func(f_min, objv, s);
            otherwise
                error('Undefined surrogate methods!')
        end
    end
end
