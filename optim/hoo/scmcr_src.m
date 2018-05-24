function [x,fval,meta] = scmcr_src(fun,x0,options)
% Method for unconstrained optimization via a cubic-regularization approach
% using a separable cubic model. Based on [Cubic-regularization counterpart
% of a variable-norm trust-region method for unconstrained minimization.
% Martínez, Raydan. 2015].
%
% Input:
%   fun     : function pointer to a function [fval,g,H] = fun(x) where fval
%             is the function value, g the gradient, and H the hessian
%             matrix at x
%   x0      : start parameters for local search
%   options : algorithm options
%
% Output:
%   x       : parameters for best found function value
%   fval    : best found function value, fun(x)
%   meta    : additional run information

% initialize values
x = x0(:);
n = size(x,1);

% paramters
if nargin < 3, options = struct(); end
options = get_options(n,options);
meta.options = options;

% function value and derivatives at starting point
[fval,g,H] = fun(x);
% initialize meta
meta.g = g;
meta.H = H;
meta.iter = 0;
meta.funEvals = 1;
% check if function differentiable at starting point
if ~isfinite(fval) || ~all(isfinite(g)) || ~all(all(isfinite(H)))
    meta.exitflag = -1;
    return;
end

% initialize further starting point variables

% [Q,T] = schur(A) where A = Q*T*Q', Q unitary, T (Schur form)
% upper triangular. Can be computed e.g. via QR. For symmetric matrices,
% this is the spectral decomposition with T = D diagonal.
[Q,D]       = schur(H);
b           = Q'*g;
gnorm       = norm(g,2);
sigma       = options.sigma0;
rho         = options.rho0;
roundoff    = sqrt(eps); % unit roundoff
y           = zeros(n,1); % space for solution of rotated problem

% main loop
while gnorm > options.tol ...
        && meta.iter < options.maxIter ...
        && meta.funEvals < options.maxFunEvals
    
    % increment iteration counter
    meta.iter = meta.iter + 1;
    
    % output
    if options.verbosity ~= 0
        if mod(meta.iter,10) == 1
            fprintf('it.\tev.\tfval\n');
        end
        fprintf('%d\t%d\t%.6e\n',meta.iter,meta.funEvals,fval);
    end
    
    % solve cubic-regularization subproblem
    for j=1:n
        c0=0;c1=b(j);c2=D(j,j)/2;c3=rho(j)/6;c4=sigma/6;
        y(j) = minPoly(c0,c1,c2,c3,c4,-options.Delta,options.Delta);
    end
    
    % compute step
    s = Q*y;
    x_new = x + s;
    fval_new = fun(x_new);
    
    % evaluate step
    
    if ~isfinite(fval_new) ...
            || fval_new > fval - options.alpha * sum(abs(y).^3)
        % no sufficient decrease: increase sigma
        sigma = max([options.sigma_small,options.eta*sigma]);
        continue;
    end
    
    % move to better x, so also compute derivatives
    [fval_new,g_new,H_new] = fun(x_new);
    meta.funEvals = meta.funEvals + 1;
    
    % check validity
    if ~isfinite(fval_new) || ~all(isfinite(g_new)) || ~all(all(isfinite(H_new)))
        % treat like failed step
        sigma = max([options.sigma_small,options.eta*sigma]);
        continue;
    end
    
    % update values
    [Q_new,D_new] = schur(H_new);
        
    % update rho inspired by a third-order secant equation
    den = Q_new'*s;
    for j=1:n
        den_j = den(j);
        if -roundoff < den_j && den_j <= 0
            den(j) = -roundoff;
        elseif 0 < den_j && den_j < roundoff
            den(j) = roundoff;
        end
    end
    rho = diag((D_new-Q_new'*H*Q_new))./den;
        
    % keep rho within bounds
    rho = min(max(rho,options.rho_min),options.rho_max);
        
    % update running variables
    x = x_new;
    fval = fval_new;
    g = g_new;
    H = H_new;
    Q = Q_new;
    D = D_new;
    gnorm = norm(g,2);
        
    b = Q'*g;
end

% fill up meta information
if gnorm >= tol % no convergenc
    meta.exitflag = 0;
else
    meta.exitflag = 1;
end
meta.algorithm = 'scmcr_src';
meta.g = g;
meta.H = H;

end % function


%% Helper functions


function z = minPoly(c0,c1,c2,c3,c4,DeltaNeg,DeltaPos)
% compute the minimum of the function h in [DeltaNeg,DeltaPos]
h = @(z) c0 + c1*z + c2*z^2 + c3*z^3 + c4*abs(z)^3;
zs = [DeltaNeg,DeltaPos];
if c2 == 0 && c3 == 0 && c4 == 0
    z = argmin(zs,h);
elseif c2 ~= 0 && c3 == 0 && c4 == 0
    zcrt = -c1/(2*c2);
    if DeltaNeg < zcrt && zcrt < DeltaPos
        zs = [zs zcrt];
    end
    z = argmin(zs,h);
elseif c3 ~= 0 && c4 == 0
    xi = 4*c2^2-12*c3*c1;
    if xi < 0
        z = argmin(zs,h);
    else
        zcrt1 = (sqrt(xi)-2*c2)/(6*c3);
        zcrt2 = (sqrt(xi)+2*c2)/(6*c3);
        if DeltaNeg < zcrt1 && zcrt1 < DeltaPos
            zs = [zs zcrt1];
        end
        if DeltaNeg < zcrt2 && zcrt2 < DeltaPos
            zs = [zs zcrt2];
        end
        z = argmin(zs,h);
    end
elseif c4 ~= 0
    zneg = minPoly(c0,c1,c2,c3-c4,0,DeltaNeg,0);
    zpos = minPoly(c0,c1,c2,c3+c4,0,0,DeltaPos);
    z = argmin([zneg,zpos],h);
end

end


function z = argmin(zs,fun)
% value z in zs such that fun(z) is minimal among all z in zs
n = size(zs,2);
z = zs(1);
fval = fun(z);
for j=2:n
    z_new = zs(j);
    fval_new = fun(z_new);
    if fval_new < fval
        z = z_new;
        fval = fval_new;
    end
end

end


function [options] = get_options(n,options_in)
% fill options with default values, and check validity

options = struct();

% default options
options.rho0        = 1*ones(n,1);
options.rho_max     = 1e3;
options.rho_min     = -options.rho_max;
options.Delta       = 2;
options.alpha       = 1e-4;
options.eta         = 10;
options.sigma0      = 0;
options.sigma_small = 0.1;
options.tol         = 1e-8;

% additional options
options.maxIter     = Inf;
options.maxFunEvals = Inf;
options.verbosity   = 1;

% fill from input
cell_fieldnames = fieldnames(options);
cell_fieldnames_in = fieldnames(options_in);

for jf = 1:length(cell_fieldnames_in)
    fieldname = cell_fieldnames_in{jf};
    if ~any(strcmp(cell_fieldnames,fieldname))
        error(['Options field ' fieldname ' does not exist.']);
    end
    options.(fieldname) = options_in.(fieldname);
end

end