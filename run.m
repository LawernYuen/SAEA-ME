close all
clear 
clc

format long;
format compact;

addpath('gpml-master', 'moead', 'nsga2', 'addGP');

problems       = {'zdt1','zdt2','zdt3','zdt4'};
problem_length = length(problems);
Dimension      = 50;
Popsize        = 200;
Iteration      = 50;
totalrun       = 20;

for i = 1 : problem_length
    problem = problems{i};
    fprintf('Running on %s...\n', problem);
    for d_idx = 1 : size(Dimension, 2)
        dimension = Dimension(:, d_idx);
        for pop_idx = 1 : size(Popsize, 2)
            popsize = Popsize(:, pop_idx);
            for itr_idx = 1 : size(Iteration, 2)
                iteration = Iteration(:, itr_idx);
                for j = 1 : totalrun
                    filename = strcat(problem,'D',num2str(dimension),'R',num2str(j),'.mat');
                    if isfile(filename)
                        continue;
                    else
                        sop               = testmop(problem, dimension);
                        sop.popsize       = popsize;
                        sop.iteration     = iteration;
                        [pop, objs]       = saea(sop);

                        save(filename,'objs');
                    end
                end
            end
        end
    end
end