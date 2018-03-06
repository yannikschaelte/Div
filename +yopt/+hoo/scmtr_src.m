function [x,fval,output] = scmtr_src(fun,x0,options)
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

% initialize values
x = x0(:);
n = size(x,1);

% initialize meta data
output.exitflag = -1;
fval = nan;
output.g = nan(n,1);
output.H = nan(n,n);
output.iterations = 0;
output.funEvals = 0;

% unit roundoff
epsilon = sqrt(eps);

% parameters
if nargin < 3, options = struct(); end
options = getOptions(n,options);
output.options = options;

% check feasibility of starting point
if any(x<lb) || any(x>ub)
    return;
end

% function value and derivatives at start point
[fval,g,H] = fun(x);
jIter = 1;
jFunEvals = 1;

% check if function differentiable at start point
if ~isfinite(fval) || ~all(isfinite(g)) || ~all(all(isfinite(H)))
    return;
end

% [Q,T] = schur(A) where A = Q*T*Q', Q unitary, T (Schur form)
% upper triangular. Can be computed e.g. via QR. For symmetric matrices,
% this is the spectral decomposition with T = D diagonal.
[Q,D] = schur(H);
b = Q'*g;
gnorm = norm(g,2);
snorm = inf;
absfvaldiff = inf;

% reserve space for solution of rotated problem
y = zeros(n,1);

% wrap function with barrier
bounded_fun = @(x,jIter) bound_fun(x,fun,lb,ub,barrier,jIter,maxIter);

% main loop
% TODO also add conditions for steplength and minimum (negative) eigenvalue
while gnorm > tol && snorm > tol && absfvaldiff > tol && jIter < maxIter && jFunEvals < maxFunEvals
    jIter = jIter + 1;
    
    % solve trust-region subproblem
    for j=1:n
        c0=0;c1=b(j);c2=D(j,j)/2;c3=rho(j)/6;c4=sigma/6;
        y(j) = minPoly(c0,c1,c2,c3,c4,-Delta,Delta);
    end
    
    % try step
    s = Q*y;
    x_new = x + s;
    fval_new = bounded_fun(x_new,jIter);
    snorm = norm(s,2);
    
    fval_diff = fval_new - fval; % desired: improvement < 0
    predicted_fval_diff = g(:)'*s + 0.5*s'*H*s;
    
    predictionRatio = fval_diff / predicted_fval_diff;
%     disp(num2str(predictionRatio));
    if predictionRatio > 0.9
        Delta = min([Delta*1.1,4]);
    elseif predictionRatio < 0.25
        Delta = max([Delta*0.9,1]);
    end
    
    if ~isfinite(fval_new) || fval_diff > - alpha * sum(abs(y).^3)
        % no success: increase sigma
        sigma = max([sigma_small,sigma_factor*sigma]);
        % if somehow better, adapt
        if fval_new >= fval
            continue;
        end
    else
        sigma = sigma/10;
        absfvaldiff = abs(fval_diff);
    end
    
    % also compute derivatives
    [fval_new,g_new,H_new] = fun(x_new);
    jFunEvals = jFunEvals + 1; % funEvals means with derivatives
    
    % check validity
    if ~isfinite(fval_new) || ~all(isfinite(g_new)) || ~all(all(isfinite(H_new)))
        sigma = max([sigma_small,sigma_factor*sigma]);
    else
        % update values
        [Q_new,D_new] = schur(H_new);
        
        % update rho inspired by a third-order secant equation
        denominator = Q_new'*s;
        for j=1:n
            dj = denominator(j);
            if -epsilon < dj && dj <= 0
                denominator(j) = -epsilon;
            elseif 0 < dj && dj < epsilon
                denominator(j) = epsilon;
            end
        end
        rho = diag((D_new-Q_new'*H*Q_new))./denominator;
        
        % keep rho within bounds
        rho = min(max(rho,rhomin),rhomax);
        
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
    if mod(jIter,20) == 0
        fprintf('Iter.\tfval\tDelta\tsnorm\tpredictionRatio--------\n');
    end
    fprintf('%d\t%.15f\t%.15f\t%.15f\t%.15f\n',jIter,fval,Delta,snorm,predictionRatio);
end


% meta information
if gnorm >= tol
    % did not converge
    output.exitflag = 0;
else
    output.exitflag = 1;
end
output.algorithm = 'rsc';
output.iterations = jIter;
output.funEvals = jFunEvals;
output.g = g;
output.H = H;

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
options.delta0      = 2;
options.delta_min   = 0.1;
options.delta_max   = 1e5;
options.eta_v       = 0.9;
options.eta_s       = 0.1;
options.gamma       = 2;
options.gamma_d     = 0.5;
options.tol         = 1e-8;

% additional options
options.maxFunEvals = 1000;

% fill from input
cell_fieldnames = fieldnames(options);
cell_fieldnames_in = cell_fieldnames(options_in);

for jf=1:length(cell_fieldnames_in)
    fieldname = cell_fieldnames_in{jf};
    if ~any(strcmp(cell_fieldnames,fieldname))
        error(['Options field ' fieldname ' does not exist.']);
    end
    options.(fieldname) = options_in.(fieldname);
end

end