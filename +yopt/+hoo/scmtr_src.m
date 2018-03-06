function [x,fval,meta] = scmtr_src(fun,x0,options)
% Method for unconstrained optimization via a trust-region approach using a
% separable cubic model. Based on [Separable cubic modeling and a
% trust-region strategy for unconstrained minimization with impact in lobal
% optimization. Martínez, Raydan. 2015].
%
% Input:
%   fun     : function pointer to a function [fval,g,H] = fun(x) where fval
%             is the function value, g the gradient, and H the hessian
%             matrix at x
%   x0      : start parameters for local search
%   options : algorithm options

% unit roundoff
roundoff = sqrt(eps);

% initialize values
x = x0(:);
n = size(x,1);

% parameters
if nargin < 3, options = struct(); end
options = getOptions(n,options);
meta.options = options;

% function value and derivatives at start point
[fval,g,H] = fun(x);
% initialize output
meta.g = g;
meta.H = H;
meta.iter = 0;
meta.funEvals = 1;
% check if function differentiable at start point
if ~isfinite(fval) || ~all(isfinite(g)) || ~all(all(isfinite(H)))
    meta.exitflag = -1;
    return;
end

% [Q,T] = schur(A) where A = Q*T*Q', Q unitary, T (Schur form)
% upper triangular. Can be computed e.g. via QR. For symmetric matrices,
% this is the spectral decomposition with T = D diagonal.
[Q,D] = schur(H);
b = Q'*g;
gnorm = norm(g,2);
delta = options.delta0;
rho = options.rho0;
% reserve space for solution of rotated problem
y = zeros(n,1);

% main loop
while gnorm > options.tol ...
        && meta.iter < options.maxIter ...
        && meta.funEvals < options.maxFunEvals
    
    % increment iteration counter
    meta.iter = meta.iter + 1;
    
    % output
    if options.verbosity ~= 0
        if mod(meta.iter,20) == 0
            fprintf('iter\tfunEvals\tfval--------\n');
        end
        fprintf('%d\t%d\t%.6e\n',meta.iter,meta.funEvals,fval);
    end
    
    % solve trust-region subproblem
    for j=1:n
        c0=0;c1=b(j);c2=D(j,j)/2;c3=rho(j)/6;
        y(j) = minPoly(c0,c1,c2,c3,0,-delta,delta);
    end
    
    % compute step
    s = Q*y;
    fval_new = fun(x+s);
    
    % evaluate step
    fval_diff = fval - fval_new;
    q = fval + g'*s + 0.5*s'*H*s;
    pred_diff = fval - q;
    ratio = fval_diff / pred_diff;
    
    % evaluate ratio
    if ratio >= options.eta_v
        x_new = x + s;
        delta_new = options.gamma * delta;
    elseif ratio >= options.eta_s
        x_new = x + s;
        delta_new = delta;
    else
        delta = options.gamma_d * delta;
        continue;
    end
    
    % x has changed, so also compute derivatives
    [fval_new,g_new,H_new] = fun(x);
    meta.funEvals = meta.funEvals + 1;
    
    % check validity
    if ~isfinite(fval_new) || ~all(isfinite(g_new)) || ~all(all(isfinite(H_new)))
        % treat like failed step
        delta = options.gamma_d * delta;
    else
        % update values
        [Q_new,D_new] = schur(H_new);
        
        % update rho inspired by a third-order secant equation
        denominator = Q_new'*s;
        for j=1:n
            denominator_j = denominator(j);
            if -roundoff < denominator_j && denominator_j <= 0
                denominator(j) = -roundoff;
            elseif 0 < denominator_j && denominator_j < roundoff
                denominator(j) = roundoff;
            end
        end
        rho = diag((D_new-Q_new'*H*Q_new))./denominator;
        
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
        delta = delta_new;
        
        b = Q'*g;
    end
    
end

% meta information
if gnorm >= tol
    % did not converge
    meta.exitflag = 0;
else
    meta.exitflag = 1;
end
meta.algorithm = 'scmtr_src';
meta.g = g;
meta.H = H;

end


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


function [options] = getOptions(n,options_in)

options = struct();

% default options
options.rho0        = 1*ones(n,1);
options.rho_max     = 1e3;
options.rho_min     = -options.rho_max;
options.delta0      = 2;
options.delta_min   = 0.1;
options.delta_max   = 1e5;
options.eta_v       = 0.9;
options.eta_s       = 0.1;
options.gamma       = 2;
options.gamma_d     = 0.5;
options.tol         = 1e-8;

% additional options
options.maxIter     = Inf;
options.maxFunEvals = inf;
options.verbosity   = 1;

% fill from input
cell_fieldnames = fieldnames(options);
cell_fieldnames_in = fieldnames(options_in);

for jf=1:length(cell_fieldnames_in)
    fieldname = cell_fieldnames_in{jf};
    if ~any(strcmp(cell_fieldnames,fieldname))
        error(['Options field ' fieldname ' does not exist.']);
    end
    options.(fieldname) = options_in.(fieldname);
end

end