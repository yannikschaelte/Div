function [] = test_pesto_jakstat(optimizer)

rng(0);

exdir=fileparts(which('test_jakstat.m'));
addpath('../../NOODLES');
addpath('../testmodels/jakstat');
load('data_jakstat','D');
fun = @(x) nllh_jakstat(x,D);
parameters = get_parameters_jakstat();

options = PestoOptions();
options.n_starts = 20;
options.mode = 'text';
if contains(optimizer,'noodles')
    options.localOptimizer = 'noodles';
else
    options.localOptimizer = optimizer;
end
options.localOptimizerOptions = local_options(optimizer,parameters);
options.obj_type = 'negative log-posterior';

used_time = cputime;
parameters_res = getMultiStarts(parameters,fun,options);
used_time = cputime - used_time;

save(fullfile(exdir, [ 'test_pesto_jakstat_' optimizer '.mat']));

end % function

function [local_options] = local_options(optimizer,parameters)

lb = parameters.min;
ub = parameters.max;

tolFun  = 1e-10;
tolX    = 1e-12;
maxFunEvals = 2000;
maxIter     = maxFunEvals;

switch optimizer
    case 'fmincon'
%         local_options = optimoptions(@fmincon,'Algorithm','interior-point');
%         local_options.Display = 'iter';
%         local_options.MaxFunctionEvaluations = maxFunEvals;
%         local_options.SpecifyObjectiveGradient = true;
%         hessian = @(x,lambda) hessianFcn(fun,x,lambda);
%         local_options.HessianFcn = hessian;
        local_options = optimoptions(@fmincon,'Algorithm','trust-region-reflective');
        local_options.Display = 'iter';
        local_options.MaxFunctionEvaluations = maxFunEvals;
        local_options.MaxIterations = maxIter;
        local_options.HessianFcn = 'objective'; % as 3rd output of fun
%         options.SubproblemAlgorithm = 'cg'; % 'cg' 'factorization'
        local_options.SpecifyObjectiveGradient = true;
    case 'scmtr_src'
        local_options = struct();
        local_options.maxFunEvals = maxFunEvals;
    case 'scmcr_src'
        local_options = struct();
        local_options.maxFunEvals = maxFunEvals;
    case 'scmcr'
        local_options = struct();
        local_options.Lb = lb;
        local_options.Ub = ub;
        local_options.MaxFunEvals = maxFunEvals;
        local_options.MaxIter = maxIter;
    case {'noodles-str','noodles-scr'}
        local_options = struct();
        local_options.lb = lb;
        local_options.ub = ub;
        local_options.derivative_fun = @noodles.NoodleProblem.objective;
        local_options.feval_max = maxFunEvals;
        switch optimizer
            case 'noodles_str'
                local_options.subproblem = noodles.SubproblemStr();
            case 'noodles_scr'
                local_options.subproblem = noodles.SubproblemScr();
        end
end

end