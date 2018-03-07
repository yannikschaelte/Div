function [] = test_pesto_jakstat(optimizer)

rng(0);

exdir=fileparts(which('test_jakstat.m'));
addpath('..');
addpath('../testmodels/jakstat');
load('data_jakstat','D');
fun = @(x) nllh_jakstat(x,D);
parameters = get_parameters_jakstat();

options = PestoOptions();
options.n_starts = 20;
options.mode = 'text';
options.localOptimizer = optimizer;
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
        local_options = optimset;
        local_options.MaxFunEvals = maxFunEvals;
        local_options.MaxIter = maxIter;
        local_options.Display = 'iter';
        local_options.Algorithm = 'interior-point';
        local_options.TolFun = tolFun;
        local_options.TolX = tolX;
        local_options.GradObj = 'on';
        local_options.PrecondBandWidth = Inf;
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
end

end