function problem = testproblems(testname, dimension)
    
    problem = struct('name', [], 'pd', [], 'domain', [], 'func', []);
    switch lower(testname)
        case 'ellipsoid'
            problem = ellipsoid(problem, dimension);
        case 'sphere'
            problem = sphere(problem, dimension);
        case 'step'
            problem = step(problem, dimension);
        case 'ackley'
            problem = ackley(problem, dimension);
        case 'rosenbrock'
            problem = rosenbrock(problem, dimension);
        case 'rastrigin'
            problem = rastrigin(problem, dimension);
        case 'branin'
            problem = branin(problem, dimension);
        case 'griewank'
            problem = griewank(problem, dimension);
        case 'schwefel'
            problem = schwefel(problem, dimension);
        case 'weierstrass'
            problem = weierstrass(problem, dimension);
        case 'schaffers_f7'
            problem = schaffers_f7(problem, dimension);
        case 'six_hump'
            problem = six_hump(problem, dimension);
        case 'f1'
            problem = f1(problem, dimension);
        case 'f4'
            problem = f4(problem, dimension);
        case 'f5'
            problem = f5(problem, dimension);
        case 'f8'
            problem = f8(problem, dimension);
        case 'f9'
            problem = f9(problem, dimension);
        case 'f12'
            problem = f12(problem, dimension);
        case 'f13'
            problem = f13(problem, dimension);
        case 'elliptic'
            problem = elliptic(problem, dimension);
        otherwise
            error('Undefined test problem name');
    end
        
end

%% Ellipsoid function generator
function p = ellipsoid(p, dim)
    
    p.name   = 'Ellipsoid';
    p.pd     = dim;
    p.domain = [-5.12 * ones(1, dim); 5.12 * ones(1, dim)];
    p.func   = @evaluate;
    
    % Ellipsoid evaluation function
    function y = evaluate(x)
        idx        = 1 : dim;
%         idx_matrix = repmat(idx, size(x, 1), 1);
%         y          = sum(idx_matrix .* x.^2, 2);
        idx = sort(idx, 'descend');
        y = sum(idx .* (x.^2), 2);
    end
end

%% Sphere
function p = sphere(p, dim)

    p.name   = 'Sphere';
    p.pd     = dim;
    p.domain = [-5.12 * ones(1, dim); 5.12 * ones(1, dim)];
    p.func   = @evaluate;
    
    % Sphere evaluation function
    function y = evaluate(x)
        y = sum(x.^2, 2);
    end
end

%% Step
function p = step(p, dim)

    p.name   = 'Step';
    p.pd     = dim;
    p.domain = [-5.12 * ones(1, dim); 5.12 * ones(1, dim)];
    p.func   = @evaluate;
    
    % Step evaluation function
    function y = evaluate(x)
        x_flr = floor(x);
        y     = sum(x_flr.^2, 2);
    end
end

%% Ackley
function p = ackley(p, dim)

    p.name   = 'Ackley';
    p.pd     = dim;
    p.domain = [-5 * ones(1, dim); 5 * ones(1, dim)];
    p.func   = @evaluate;
    
    % Ackley evaluation function
    function y = evaluate(x)
        a = 20 * exp(-0.2 * sqrt(1 / dim * sum(x.^2, 2)));
        b = exp(1 / dim * sum(cos(2 * pi * x), 2));
        y = 20 + exp(1) - a - b;
    end   
end

%% Rosenbrock
function p = rosenbrock(p, dim)

    p.name   = 'Rosenbrock';
    p.pd     = dim;
    p.domain = [-2.048 * ones(1, dim); 2.048 * ones(1, dim)];
    p.func   = @evaluate;
    
    % Rosenbrock evaluation function
    function y = evaluate(x)
        xi  = x(:, 1 : dim-1);
        xi1 = x(:, 2 : dim);
        yi  = 100 * (xi1 - xi.^2) .^2 + (xi - 1) .^2;
        y   = sum(yi, 2);
    end   
end

%% Rastrigin
function p = rastrigin(p, dim)

    p.name   = 'Rastrigin';
    p.pd     = dim;
    p.domain = [-5.12 * ones(1, dim); 5.12 * ones(1, dim)];
    p.func   = @evaluate;
    
    % Rastrigin evaluation function
    function y = evaluate(x)
        yi = x .^2 - 10 * cos(2 * pi * x);
        y  = 10 * dim + sum(yi, 2);
    end   
end

%% Griewank
function p = griewank(p, dim)

    p.name   = 'Griewank';
    p.pd     = dim;
    p.domain = [-50 * ones(1, dim); 50 * ones(1, dim)];
    p.func   = @evaluate;
    
    % Griewank evaluation function
    function y = evaluate(x)
        y1  = sum(x.^2, 2) / 4000;
        idx = 1 : dim;
        y2  = prod(cos(x./sqrt(idx)), 2);
        y   = y1 - y2 + 1;
    end   
end

%% Schwefel
function p = schwefel(p, dim)

    p.name   = 'Schwefel';
    p.pd     = dim;
    p.domain = [-500 * ones(1, dim); 500 * ones(1, dim)];
    p.func   = @evaluate;
    
    % Schwefel evaluation function
    function y = evaluate(x)
        y1 = x .* sin(sqrt(abs(x)));
        y  = 418.9829 * dim - sum(y1, 2);
    end   
end

%% Weierstrass
function p = weierstrass(p, dim)

    p.name   = 'Weierstrass';
    p.pd     = dim;
    p.domain = [-50 * ones(1, dim); 50 * ones(1, dim)];
    p.func   = @evaluate;
    
    % Weierstrass evaluation function
    function y = evaluate(x)
        k  = 0:20;
        a  = 0.5.^k;
        b  = 3.^k;
%         y1 = bsxfun(@times, x+0.5, reshape(b,1,1,numel(b)));
%         y2 = cos(2 * pi * y1);
%         y3 = bsxfun(@times, y2, reshape(a,1,1,[]));
%         y4 = sum(sum(y3, 3), 2);
%         y5 = a .* cos(pi * b);
%         y6 = dim * sum(y5, 2);
%         y  = y4 - y6;
        y1 = 0;
        for i = 1 : numel(a)
            y1 = y1 + a(i) * cos(2*pi*b(i)*(x+0.5));
        end
        y2 = sum(y1, 2);
        y3 = a .* cos(pi * b);
        y4 = dim * sum(y3, 2);
        y  = y2 - y4;
    end
end

%% Schaffers F7
function p = schaffers_f7(p, dim)

    p.name   = 'Schaffers_F7';
    p.pd     = dim;
    p.domain = [-30 * ones(1, dim); 30 * ones(1, dim)];
    p.func   = @evaluate;
    
    % Schaffers F7 evaluation function
    function y = evaluate(x)
        x1 = x(:, 1:end-1);
        x2 = x(:, 2:end);
        z  = sqrt(x1.^2 + x2.^2);
        y1 = sqrt(z) + sqrt(z) .* (sin(50*z.^0.2).^2);
        y2 = sum(y1, 2) / (dim-1);
        y  = y2 .^ 2;
    end
end

%% Six-Hump Camel function generator
function p = six_hump(p, dim)

    p.name   = 'Six_Hump_Camel';
    p.pd     = dim;
    p.domain = [-3.0, -2.0; 3.0, 2.0];
    p.func   = @evaluate;
    
    % Six-Hump Camel evaluation function
    function y = evaluate(x)
        x1 = x(:, 1);
        x2 = x(:, 2);
        y1 = 4 - 2.1*x1.^2 + x1.^4/3;
        y2 = (4*x2.^2 - 4) .* x2.^2;
        y  = y1.*x1.^2 + x1.*x2 + y2;
    end
end

%% Branin function generator
function p = branin(p, dim)

    p.name   = 'Branin';
    p.pd     = dim;
    p.domain = [-5.0, 0.0; 10.0, 15.0];
    p.func   = @evaluate;
    
    % Branin evaluation function
    function y = evaluate(x)
        a = 1; b = 5.1 / (4 * pi * pi); c = 5 / pi; d = 6; h = 10;
        ff = 1 / (8 * pi);
        x1 = x(:, 1);
        x2 = x(:, 2);
        y  = a .* (x2 - b .* x1.^2 + c .* x1 - d).^2 + ...
            h .* (1 - ff) .* cos(x1) + h;        
    end
end

%% Large-scale problems
%------------------------------------------------------------------------------
% Elliptic Function
%------------------------------------------------------------------------------
function fit = elliptic(x)

    D = size(x, 2);
    condition = 1e+6;
    coefficients = condition .^ linspace(0, 1, D); 
    fit = T_irreg(x).^2 * coefficients'; 
end

%------------------------------------------------------------------------------
% Rastrigin's Function
%------------------------------------------------------------------------------
function fit = rastrigins(x)
    D = size(x, 2);
    A = 10;
    x = T_diag(T_asy(T_irreg(x), 0.2), 10);
    fit = A*(D - sum(cos(2*pi*x), 2)) + sum(x.^2, 2);
end

%------------------------------------------------------------------------------
% Schwefel's Problem 1.2
%------------------------------------------------------------------------------
function fit = schwefels(x)
    D = size(x, 2);
    x = T_asy(T_irreg(x), 0.2);
    fit = 0;
    for i = 1:D
        fit = fit + sum(x(:,1:i),2).^2;
    end
end

%------------------------------------------------------------------------------
% Rosenbrock's Function
%------------------------------------------------------------------------------
function fit = rosenbrocks(x)
    D = size(x, 2);
    fit = sum(100.*(x(:,1:D-1).^2-x(:,2:D)).^2+(x(:,1:D-1)-1).^2, 2);
end

%------------------------------------------------------------------------------
% f1: Shifted Elliptic Function
% D = 100
%------------------------------------------------------------------------------
function p = f1(p, dim)

    assert(dim==100, 'The dimension should be 100 in this test function.');
    p.name   = 'f1';
    p.pd     = dim;
    p.domain = [-10 * ones(1, dim); 10 * ones(1, dim)];
    p.func   = @evaluate;

    function y = evaluate(x)
        idx = checkBounds(x, p.domain);
        y = elliptic(x);
        y(idx) = NaN;
        if ~isempty(idx)
            warning 'Some of the solutions are violating boundary constraints.';
        end
    end
end

%------------------------------------------------------------------------------
% f4: 7-nonseparable, 1-separable Shifted and Rotated Elliptic Function
% D = 100
%------------------------------------------------------------------------------
function p = f4(p, dim)

    assert(dim==100, 'The dimension should be 100 in this test function.');
    p.name   = 'f4';
    p.pd     = dim;
    p.domain = [-10 * ones(1, dim); 10 * ones(1, dim)];
    p.func   = @evaluate;
    p.s      = [10;5;5;15;10;5;5;45];
    if ~isfield(p,'permut') || isempty(p.permut)
        p.permut = randperm(dim);
    end
    if ~isfield(p,'ro5') || isempty(p.ro5)
       p.ro5    = orth(rand(5,5));
    end
    if ~isfield(p,'ro10') || isempty(p.ro10)
       p.ro10   = orth(rand(10,10));
    end
    if ~isfield(p,'ro15') || isempty(p.ro15)
      p.ro15   = orth(rand(15,15));
    end
    if ~isfield(p,'w') || isempty(p.w)
      p.w      = 10.^(3*randn(1,7));
    end

    function y = evaluate(x)
        idx = checkBounds(x, p.domain);
        y = 0;
        ldim = 1;
        for i=1:length(p.s)-1
            if (p.s(i) == 5)
                f = elliptic(x(:, p.permut(ldim:ldim+p.s(i)-1))*p.ro5);
                ldim = ldim + p.s(i);
            elseif (p.s(i) == 10)
                f = elliptic(x(:, p.permut(ldim:ldim+p.s(i)-1))*p.ro10);
                ldim = ldim + p.s(i);
            elseif (p.s(i) == 15)
                f = elliptic(x(:, p.permut(ldim:ldim+p.s(i)-1))*p.ro15);
                ldim = ldim + p.s(i);
            end
            y = y + p.w(i)*f;
        end
        y = y + elliptic(x(:, p.permut(ldim:end)));
        y(idx) = NaN;
        if ~isempty(idx)
            warning 'Some of the solutions are violating boundary constraints.';
        end        
    end
end

%------------------------------------------------------------------------------
% f5: 7-nonseparable, 1-separable Shifted and Rotated Rastrigin’s Function
% D = 100
%------------------------------------------------------------------------------
function p = f5(p, dim)
    
    assert(dim==100, 'The dimension should be 100 in this test function.');
    p.name   = 'f5';
    p.pd     = dim;
    p.domain = [-5 * ones(1, dim); 5 * ones(1, dim)];
    p.func   = @evaluate;
    p.s      = [10;5;5;15;10;5;5;45];
    if ~isfield(p,'permut') || isempty(p.permut)
        p.permut = randperm(dim);
    end
    if ~isfield(p,'ro5') || isempty(p.ro5)
       p.ro5    = orth(rand(5,5));
    end
    if ~isfield(p,'ro10') || isempty(p.ro10)
       p.ro10   = orth(rand(10,10));
    end
    if ~isfield(p,'ro15') || isempty(p.ro15)
      p.ro15   = orth(rand(15,15));
    end
    if ~isfield(p,'w') || isempty(p.w)
      p.w      = 10.^(3*randn(1,7));
    end
    
    function y = evaluate(x)
        idx = checkBounds(x, p.domain);
        y = 0;
        ldim = 1;
        for i=1:length(p.s)-1
            if (p.s(i) == 5)
                f = rastrigins(x(:, p.permut(ldim:ldim+p.s(i)-1))*p.ro5);
                ldim = ldim + p.s(i);
            elseif (p.s(i) == 10)
                f = rastrigins(x(:, p.permut(ldim:ldim+p.s(i)-1))*p.ro10);
                ldim = ldim + p.s(i);
            elseif (p.s(i) == 15)
                f = rastrigins(x(:, p.permut(ldim:ldim+p.s(i)-1))*p.ro15);
                ldim = ldim + p.s(i);
            end
            y = y + p.w(i)*f;
        end
        y = y + rastrigins(x(:, p.permut(ldim:end)));
        y(idx) = NaN;
        if ~isempty(idx)
            warning 'Some of the solutions are violating boundary constraints.';
        end
    end
end

%------------------------------------------------------------------------------
% f8: 11-nonseparable Shifted and Rotated Elliptic Function
% D = 100
%------------------------------------------------------------------------------
function p = f8(p, dim)
    
    assert(dim==100, 'The dimension should be 100 in this test function.');
    p.name   = 'f8';
    p.pd     = dim;
    p.domain = [-10 * ones(1, dim); 10 * ones(1, dim)];
    p.func   = @evaluate;
    p.s      = [10;10;5;5;15;15;5;5;10;5;15];
    if ~isfield(p,'permut') || isempty(p.permut)
        p.permut = randperm(dim);
    end
    if ~isfield(p,'ro5') || isempty(p.ro5)
       p.ro5    = orth(rand(5,5));
    end
    if ~isfield(p,'ro10') || isempty(p.ro10)
       p.ro10   = orth(rand(10,10));
    end
    if ~isfield(p,'ro15') || isempty(p.ro15)
      p.ro15   = orth(rand(15,15));
    end
    if ~isfield(p,'w') || isempty(p.w)
      p.w      = 10.^(3*randn(1,11));
    end
    
    function y = evaluate(x)
        idx = checkBounds(x, p.domain);
        y = 0;
        ldim = 1;
        for i=1:length(p.s)
            if (p.s(i) == 5)
                f = elliptic(x(:, p.permut(ldim:ldim+p.s(i)-1))*p.ro5);
                ldim = ldim + p.s(i);
            elseif (p.s(i) == 10)
                f = elliptic(x(:, p.permut(ldim:ldim+p.s(i)-1))*p.ro10);
                ldim = ldim + p.s(i);
            elseif (p.s(i) == 15)
                f = elliptic(x(:, p.permut(ldim:ldim+p.s(i)-1))*p.ro15);
                ldim = ldim + p.s(i);
            end
            y = y + p.w(i)*f;
        end
        y(idx) = NaN;
        if ~isempty(idx)
            warning 'Some of the solutions are violating boundary constraints.';
        end
    end
end

%------------------------------------------------------------------------------
% f9: 11-nonseparable Shifted and Rotated Rastrigin’s Function
% D = 100
%------------------------------------------------------------------------------
function p = f9(p, dim)
    
    assert(dim==100, 'The dimension should be 100 in this test function.');
    p.name   = 'f9';
    p.pd     = dim;
    p.domain = [-5 * ones(1, dim); 5 * ones(1, dim)];
    p.func   = @evaluate;
    p.s      = [10;5;5;15;10;5;5;45];
    if ~isfield(p,'permut') || isempty(p.permut)
        p.permut = randperm(dim);
    end
    if ~isfield(p,'ro5') || isempty(p.ro5)
       p.ro5    = orth(rand(5,5));
    end
    if ~isfield(p,'ro10') || isempty(p.ro10)
       p.ro10   = orth(rand(10,10));
    end
    if ~isfield(p,'ro15') || isempty(p.ro15)
      p.ro15   = orth(rand(15,15));
    end
    if ~isfield(p,'w') || isempty(p.w)
      p.w      = 10.^(3*randn(1,11));
    end
    
    function y = evaluate(x)
        idx = checkBounds(x, p.domain);
        y = 0;
        ldim = 1;
        for i=1:length(p.s)
            if (p.s(i) == 5)
                f = rastrigins(x(:, p.permut(ldim:ldim+p.s(i)-1))*p.ro5);
                ldim = ldim + p.s(i);
            elseif (p.s(i) == 10)
                f = rastrigins(x(:, p.permut(ldim:ldim+p.s(i)-1))*p.ro10);
                ldim = ldim + p.s(i);
            elseif (p.s(i) == 15)
                f = rastrigins(x(:, p.permut(ldim:ldim+p.s(i)-1))*p.ro15);
                ldim = ldim + p.s(i);
            end
            y = y + p.w(i)*f;
        end
        y(idx) = NaN;
        if ~isempty(idx)
            warning 'Some of the solutions are violating boundary constraints.';
        end
    end
end

%------------------------------------------------------------------------------
% f12: Shifted Rosenbrock’s Function
% D = 100
%------------------------------------------------------------------------------
function p = f12(p, dim)

    assert(dim==100, 'The dimension should be 100 in this test function.');
    p.name   = 'f12';
    p.pd     = dim;
    p.domain = [-10 * ones(1, dim); 10 * ones(1, dim)];
    p.func   = @evaluate;

    function y = evaluate(x)
        idx = checkBounds(x, p.domain);
        y = rosenbrocks(x);
        y(idx) = NaN;
        if ~isempty(idx)
            warning 'Some of the solutions are violating boundary constraints.';
        end
    end
end

%------------------------------------------------------------------------------
% f13: Shifted Schwefel’s Function with Conforming Overlapping Subcomponents
% D = 100
%------------------------------------------------------------------------------
function p = f13(p, dim)

    assert(dim==100, 'The dimension should be 100 in this test function.');
    p.name   = 'f13';
    p.pd     = dim;
    p.domain = [-10 * ones(1, dim); 10 * ones(1, dim)];
    p.func   = @evaluate;
    p.s      = [10;5;5;15;10;5;5;45];
    p.m      = 5;
    if ~isfield(p,'permut') || isempty(p.permut)
        p.permut = randperm(dim);
    end
    if ~isfield(p,'ro5') || isempty(p.ro5)
       p.ro5    = orth(rand(5,5));
    end
    if ~isfield(p,'ro10') || isempty(p.ro10)
       p.ro10   = orth(rand(10,10));
    end
    if ~isfield(p,'ro15') || isempty(p.ro15)
      p.ro15   = orth(rand(15,15));
    end
    if ~isfield(p,'w') || isempty(p.w)
      p.w      = 10.^(3*randn(1,11));
    end

    function y = evaluate(x)
        idx = checkBounds(x, p.domain);
        y = 0;
        ldim = 1;
        c = cumsum(p.s);
        for i=1:length(p.s)
            if i ~= 1
                ldim = c(i-1) - ((i-1)*p.m) + 1;
            end
            udim = c(i) - ((i-1)*p.m);
            if (p.s(i) == 5)
                f = schwefels(x(:, p.permut(ldim:udim))*p.ro5);
            elseif (p.s(i) == 10)
                f = schwefels(x(:, p.permut(ldim:udim))*p.ro10);
            elseif (p.s(i) == 15)
                f = schwefels(x(:, p.permut(ldim:udim))*p.ro15);
            end
            y = y + p.w(i)*f;
        end
        y(idx) = NaN;
        if ~isempty(idx)
            warning 'Some of the solutions are violating boundary constraints.';
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------------------------------------------------
% This transformation function is used to break the symmetry of symmetric 
% functions.
%------------------------------------------------------------------------------
function g = T_asy(f, beta)
    [popsize, D] = size(f);
    g = f;
    temp = repmat(beta * linspace(0, 1, D), popsize, 1); 
    ind = f > 0;
    g(ind) = f(ind).^ (1 + temp(ind) .* sqrt(f(ind)));  
end


%------------------------------------------------------------------------------
% This transformation is used to create the ill-conditioning effect.
%------------------------------------------------------------------------------
function g = T_diag(f, alpha)
    [popsize, D] = size(f);
    scales = repmat(sqrt(alpha) .^ linspace(0, 1, D), popsize, 1); 
    g = scales .* f;
end


%------------------------------------------------------------------------------
% This transformation is used to create smooth local irregularities.
%------------------------------------------------------------------------------
function g = T_irreg(f)
   a = 0.1;
   g = f; 
   idx = (f > 0);
   g(idx) = log(f(idx))/a;
   g(idx) = exp(g(idx) + 0.49*(sin(g(idx)) + sin(0.79*g(idx)))).^a;
   idx = (f < 0);
   g(idx) = log(-f(idx))/a;
   g(idx) = -exp(g(idx) + 0.49*(sin(0.55*g(idx)) + sin(0.31*g(idx)))).^a;
end


%------------------------------------------------------------------------------
% This function tests a given decision vector against the boundaries of a function.
%------------------------------------------------------------------------------
function indices = checkBounds(x, domain)
    indices = find(sum(x > domain(2,:) | x < domain(1,:)) > 0);
end

