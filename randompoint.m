function parent_pop = randompoint(problem, popsize, type)
    % This function initializes the parent population
    %
    % problem    : The structure of the current problem
    % popsize    : The population size of the paretn population
    % type       : The sampling method
    %
    % parent_pop : The initialized parent population

    if strcmp(type, 'random')
        randarray  = rand(popsize, problem.pd);
        lowend     = problem.domain(1, :);
        span       = problem.domain(2, :) - lowend;
        parent_pop = randarray .* span + lowend;
    elseif strcmp(type, 'lhs')
        lowend     = problem.domain(1, :);
        span       = problem.domain(2, :) - lowend;
        parent_pop = lowend + span .* lhsdesign(popsize, problem.pd, 'criterion', 'maximin', 'iterations', 1000);
    else
        error('Undefined sampling method!')
    end
    clear lowend span
end
