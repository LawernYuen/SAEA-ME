function acq = acquisition_func(f_min, m, s)

    global params;

    switch params.acq_type
        case 'EI'
            s    = sqrt(s) + 1e-6;
            diff = f_min - m - 1e-6;
            norm = diff ./ s;
            EI   = diff .* normcdf(norm) + s .* normpdf(norm);
            acq  = -EI;
            
            clear diff EI norm m s
        case 'PI'
            s    = sqrt(s);
            diff = f_min - m;
            norm = diff ./ s;
            PI   = normcdf(norm);
            acq  = -PI;
        case 'UCB'
            s    = sqrt(s);
            beta = 0.5;
            UCB  = m + beta * s;
            acq  = -UCB;
        case 'LCB'
            s    = sqrt(s);
            w    = 2;
            LCB  = m - w * s;
            acq  = -LCB;
        otherwise
            error('Undefined acquisition function!') 
    end

end

function y = gaussian_PDF(x)
    y = 1 / sqrt(2 * pi) * exp(-x.^2 / 2);
end

function y = gaussian_CDF(x)
    y = 0.5 * (1 + erf(x / sqrt(2)));
end