function [] = test_rme(startpoint_method, random_index)

addpath('../../testmodels/rafmekerk');

rng('default');
rng(random_index);

totalMaxFunEvals = 10000;
nStarts = 10;
maxFunEvals = round( totalMaxFunEvals / nStarts );

solver = 'fmincon';

load('data_rafmekerk_noreps.mat', 'D');
objfun = @(x) llh_rafmekerk_standard(x, D);

[parameters, ~] = getParametersAndOptions_rafmekerk('standard');
parameters.guess = [];

if strcmp(startpoint_method, 'latin hypercube')
    ss_maxFunEvals = 0;
else
    ss_maxFunEvals = max([2 * nStarts, min([nStarts * 50, round(totalMaxFunEvals / nStarts)])]);
end    

options = PestoOptions();
options.obj_type = 'log-posterior';
options.proposal = startpoint_method;
options.ss_maxFunEvals = ss_maxFunEvals;
options.n_starts = nStarts;
options.objOutNumber = 2;
options.mode = 'text';
options.localOptimizer = solver;

lOptions = optimoptions(@fmincon);
lOptions.MaxFunctionEvaluations = round( (totalMaxFunEvals - ss_maxFunEvals) / nStarts );
lOptions.MaxIterations = inf;
lOptions.Display = 'off';
lOptions.GradObj = 'on';
options.localOptimizerOptions = lOptions;

disp(['ss_maxFunEvals: ' num2str(options.ss_maxFunEvals)]);
disp(['local maxFunEvals: ' num2str(lOptions.MaxFunctionEvaluations)]);

starttime = tic;
parameters_res = getMultiStarts(parameters, objfun, options);
time = toc(starttime);
parameters_res.time = time;

save(['res/test_rme_' strrep(startpoint_method, ' ', '-') '_' num2str(random_index) '_' num2str(maxFunEvals) '_' num2str(nStarts)]);

end
