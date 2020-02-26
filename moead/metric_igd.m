addpath('hv-1.3-src');

% load('Data_PF/UF2.mat');
PF = new_sample_obj;
PF_size  = size(PF, 1);

Reference_point = [2.0, 2.0];

for i = 1 : 1
    load(['pareto_dtlz2_', int2str(i),'.mat'])
    
    pop = pareto;    
    popsize  = size(pareto, 1);
    min_dist = zeros(PF_size, 1);
    
    for j = 1 : PF_size
        temp    = PF(j, :);
        tempMat = temp(ones(1, popsize), :);
        
        temp_dist   = (tempMat - pop) .^ 2;
        distance    = sum(temp_dist, 2);
        min_dist(j) = min(distance);
    end
    
    min_dist = sqrt(min_dist);
    igd(i) = mean(min_dist);
    
%     HV(i) = Hypervolume_MEX(pop, Reference_point);
end

mean_igd = mean(igd)
std_igd  = std(igd)

% mean_hv = mean(HV)
% std_hv  = std(HV)