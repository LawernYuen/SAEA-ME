function ind = gaussian_mutate(ind, prob, domain)
%GAUSSIAN_MUTATE Summary of this function goes here
%   Detailed explanation goes here

   parDim  = length(ind);
   lowend  = domain(1, :);
   highend = domain(2, :);
   sigma   = (highend - lowend) ./ 20;
   
   newparam = min(max(normrnd(ind, sigma), lowend), highend);
   C        = rand(parDim, 1) < prob;
   ind(C)   = newparam(C);
end