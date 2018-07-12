addpath('../testmodels/rafmekerk');

rng('default');
rng(0);

totalMaxFunEvals = 10000;
nStarts = 10;
maxFunEvals = round( totalMaxFunEvals / nStarts );

solver = 'fmincon';

load data_rafmekerk_noreps.mat
objfun = @(x) llh_rafmekerk_standard(x, D);

[parameters, ~] = getParametersAndOptions_rafmekerk('standard');
lb = parameters.min;
ub = parameters.max;
nPar = parameters.number;

options = PestoOptions();
options.obj_type = 'log-posterior';
options.proposal = 'latin hypercube';
options.n_starts = nStarts;
% parallelized duration * 2
options.ss_maxFunEvals = 0;
% options.ss_maxFunEvals = max([2 * nStarts, min([nStarts * 50, totalMaxFunEvals / nStarts])]);
disp(['ss_maxFunEvals: ' num2str(options.ss_maxFunEvals)]);
options.objOutNumber = 2;
options.mode = 'text';
options.localOptimizer = solver;

lOptions = optimoptions(@fmincon);
lOptions.MaxFunctionEvaluations = round ( (totalMaxFunEvals - options.ss_maxFunEvals) / nStarts );
disp(['local maxFunEvals: ' num2str(lOptions.MaxFunctionEvaluations)]);
lOptions.MaxIterations = inf;
lOptions.Display = 'off';
lOptions.GradObj = 'on';
options.localOptimizerOptions = lOptions;

starttime = tic;
parameters_res = getMultiStarts(parameters, objfun, options);
time = toc(starttime);
parameters_res.time = time;

save(['test_rme_' options.proposal '_' num2str(maxFunEvals) '_' num2str(nStarts)], 'parameters_res');