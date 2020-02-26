function best_ind = cmaes(problem, model_hyp, db_pop, db_objs, f_min)   % (mu/mu_w, lambda)-CMA-ES

    global nFEs params meanfunc covfunc likfunc inf;

    %% --------------------  Initialization --------------------------------  
    % User defined input parameters (need to be edited)
    xmean = rand(1, problem.pd);    % objective variables initial point
    sigma = 0.3;          % coordinate wise standard deviation (step size)
    stopfitness = 1e-10;  % stop if fitness < stopfitness (minimization)
    stopeval = 1e3;   % stop after stopeval number of function evaluations

    % Strategy parameter setting: Selection  
    lambda = 4 + floor(3 * log(problem.pd));  % population size, offspring number
    mu = lambda / 2;               % number of parents/points for recombination
    weights = log(mu + 1/2) - log(1:mu); % muXone array for weighted recombination
    mu = floor(mu);        
    weights = weights / sum(weights);     % normalize recombination weights array
    mueff = sum(weights)^2 / sum(weights.^2); % variance-effectiveness of sum w_i x_i

    % Strategy parameter setting: Adaptation
    cc = (4 + mueff/problem.pd) / (problem.pd + 4 + 2*mueff/problem.pd);  % time constant for cumulation for C
    cs = (mueff + 2) / (problem.pd + mueff + 5);  % t-const for cumulation for sigma control
    c1 = 2 / ((problem.pd + 1.3)^2 + mueff);    % learning rate for rank-one update of C
    cmu = min(1 - c1, 2 * (mueff - 2 + 1/mueff) / ((problem.pd + 2)^2 + mueff));  % and for rank-mu update
    damps = 1 + 2*max(0, sqrt((mueff-1)/(problem.pd+1))-1) + cs; % damping for sigma 
                                                      % usually close to 1
    % Initialize dynamic (internal) strategy parameters and constants
    pc = zeros(1,problem.pd); ps = zeros(1,problem.pd);   % evolution paths for C and sigma
    B = eye(problem.pd,problem.pd);                       % B defines the coordinate system
    D = ones(1,problem.pd);                      % diagonal D defines the scaling
    C = B * diag(D.^2) * B';            % covariance matrix C
    invsqrtC = B * diag(D.^-1) * B';    % C^-1/2 
    eigeneval = 0;                      % track update of B and D
    chiN=problem.pd^0.5*(1-1/(4*problem.pd)+1/(21*problem.pd^2));  % expectation of 
                                      %   ||problem.pd(0,I)|| == norm(randn(problem.pd,1)) 
    %% -------------------- Generation Loop --------------------------------
    counteval = 0;  % the next 40 lines contain the 20 lines of interesting code 
    while counteval < stopeval

        % Generate and evaluate lambda offspring
        arx = ones(lambda, problem.pd);
        for k=1:lambda
            arx(k,:) = xmean + sigma * (D .* randn(1,problem.pd)) * B; % m + sig * Normal(0,C) 
            counteval = counteval+1;
        end
        switch params.model
            case 'vanilla'
                arfitness = problem.func(arx);  % use true obj evaluation
            case 'gp'
                [newobj, news] = gp(model_hyp, inf, meanfunc, covfunc, likfunc, db_pop, db_objs, arx);
                arfitness   = acquisition_func(f_min, newobj, news);
            case 'kriging'
                [newobj, news] = predictor(arx, model_hyp);
                arfitness = acquisition_func(f_min, newobj, news);
            otherwise
                error('Undefined model type!')
        end

        % Sort by fitness and compute weighted mean into xmean
        [arfitness, arindex] = sort(arfitness); % minimization
        xold = xmean;
        xmean = weights*arx(arindex(1:mu),:);   % recombination, new mean value

        % Cumulation: Update evolution paths
        ps = (1-cs)*ps ... 
            + sqrt(cs*(2-cs)*mueff) * (xmean-xold) / sigma * invsqrtC; 
        hsig = norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN < 1.4 + 2/(problem.pd+1);
        pc = (1-cc)*pc ...
            + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma;

        % Adapt covariance matrix C
        artmp = (1/sigma) * (arx(arindex(1:mu),:)-repmat(xold,mu,1));
        C = (1-c1-cmu) * C ...                  % regard old matrix  
            + c1 * (pc*pc' ...                 % plus rank one update
                + (1-hsig) * cc*(2-cc) * C) ... % minor correction if hsig==0
                + cmu * artmp' * diag(weights) * artmp; % plus rank mu update

        % Adapt step size sigma
        sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1)); 

        % Decomposition of C into B*diag(D.^2)*B' (diagonalization)
        if counteval - eigeneval > lambda/(c1+cmu)/problem.pd/10  % to achieve O(problem.pd^2)
            eigeneval = counteval;
            C = triu(C) + triu(C,1)'; % enforce symmetry
            [B,D] = eig(C);           % eigen decomposition, B==normalized eigenvectors
            D = sqrt(diag(D))';        % D is a vector of standard deviations now
            invsqrtC = B * diag(D.^-1) * B';
        end

        fprintf('FEs = %d :: Itr = %d :: Best EI = %f :: Best Fitness = %f\n', nFEs, counteval, arfitness(1), f_min);

        % Break, if fitness is good enough or condition exceeds 1e14, better termination methods are advisable 
%         if arfitness(1) <= stopfitness || max(D) > 1e7 * min(D)
        if max(D) > 1e7 * min(D)
            break;
        end
        
    end % while, end generation loop

    % repair the solution
    arx(arindex(1), :) = max(arx(arindex(1), :), problem.domain(1, :));
    arx(arindex(1), :) = min(arx(arindex(1), :), problem.domain(2, :));
    best_ind = arx(arindex(1), :); % Return best point of last iteration.
                             % Notice that xmean is expected to be even
                             % better.
    % ---------------------------------------------------------------  
end