
function [db_pop, db_objs] = saea(problem)

    global nFEs params;
    
    % set the algorithm parameters
    [db_pop, db_objs] = init(problem);
    
    while ~terminate(nFEs)
        tr_pop  = db_pop;
        tr_objs = db_objs;
        
        % build (or re-build) the surrogate model
        model_hyp = cell(1,problem.od);
        for i = 1:problem.od
            if strcmp(params.model, 'addGP')
                [params.decomp, params.paramsGP] = ...
                    preprocessDecomposition(numel(params.seps{i}), params, struct, true);
            end
            model_hyp{i} = model_building(problem, tr_pop(:,params.seps{i}), tr_objs(:,i));
        end
        
        % Optimise acquisition function
        switch params.optimize
            case 'moead'
                [pop, popobj] = moead(problem, model_hyp, tr_pop, tr_objs);
            case 'nsga2'
                [pop, popobj] = nsga2(problem, model_hyp, tr_pop, tr_objs);
        end
        
        % Get soultions of interest by HV
        pf = problem.pf();
        objLB = [popobj(:,2), popobj(:,4)];
        obj   = [popobj(:,1), popobj(:,3)];
        objLHV = zeros(size(pop,1), 1);
        objHV  = zeros(size(pop,1), 1);
        popLHV = HV(objLB, pf);
        popHV = HV(obj, pf);
        for i = 1:size(pop,1)
            temp = [1:i-1, i+1:size(pop,1)];
            objLHV(i) = HV(objLB(temp,:), pf);
            objHV(i) = HV(obj(temp,:), pf);
        end
        objLHV = popLHV - objLHV;
        objHV = popHV - objHV;
        [~,idx1] = sort(objLHV,'descend');
        [~,idx2] = sort(objHV,'descend');
        if size(pop,1) < 10
            soi = 1:size(pop, 1);
        else
            o1 = idx1(1:10);
            o2 = idx2(1:10);
            soi = intersect(o1, o2);
            if numel(soi) < 5
                soi = [soi; idx1(1:3); idx2(1:3)];
            end
        end
        
        new_sample = pop(soi, :);
        new_sample_obj = problem.func(new_sample);
        nFEs = nFEs + size(new_sample, 1);
        db_pop  = [db_pop; new_sample];
        db_objs = [db_objs; new_sample_obj];
        
    end
        
end

%% initialisation process
function [pop, objs] = init(problem)

    global nFEs params meanfunc covfunc likfunc inf;
    
    % parameter settings
    params.model    = 'gp'; % model type: GP, Kriging
    params.db_size  = 11 * problem.pd - 1;
    params.popsize  = problem.popsize;  % population size of DE
    params.maxFEs   = 800;
    params.maxItr   = problem.iteration;  % number of generations of DE
    params.optimize = 'nsga2';
    
    switch params.model
        case 'gp'
            % parameter settings of GP
            startup;
            meanfunc = @meanConst;
            covfunc  = @covSEiso;
            likfunc  = @likGauss;
            inf      = @infGaussLik;
        case 'kriging'
            meanfunc = @regpoly0;
            covfunc  = @corrgauss;
        case 'addGP'
    end
    
    % initialise a population
    pop  = randompoint(problem, params.db_size, 'lhs');
    objs = problem.func(pop);
    nFEs = params.db_size;
    
    [seps, popAddition, objAddition] = group(problem);
    pop = [pop; popAddition];
    objs = [objs; objAddition];
    params.seps = seps;
    
end

%% model training and hyperparameter estimation
function model = model_building(problem, db_pop, db_objs)
    global meanfunc covfunc likfunc inf params;
    
    if strcmp(params.model, 'vanilla') == 1
        model = [];
    elseif strcmp(params.model, 'gp') == 1
        hyp.mean = 0.5;
        hyp.cov  = [1; 1];
        sn      = 0.1; 
        hyp.lik = log(sn);
    
        hyp   = minimize(hyp, @gp, -200, inf, meanfunc, covfunc, likfunc, db_pop, db_objs);
        model = hyp;    % this represents the GP model hyperparameters
    elseif strcmp(params.model, 'kriging') == 1
        theta = ones(1, problem.pd);
        lob   = 0.001 * ones(1, problem.pd);
        upb   = 1000 * ones(1, problem.pd);

        model = dacefit(db_pop, db_objs, meanfunc, covfunc, theta, lob, upb);
    elseif strcmp(params.model, 'addGP') == 1
        [combFuncH, funcHs] = addGPBO(problem, params.paramsGP, params.decomp, db_pop, db_objs, params.numItr);
        model = {combFuncH, funcHs};
    else
        error('Undefined surrogate model!')
    end
    
end

%% termination checking
function y = terminate(nFEs)

    global params;
    y = nFEs >= params.maxFEs;
    
end