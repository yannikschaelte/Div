function [] = test_jakstat(optimizer)

exdir=fileparts(which('test_jakstat.m'));
addpath('..');
addpath('../testmodels/jakstat');
addpath('../../Noodles');
load('data_jakstat','D');
fun = @(x) nllh_jakstat(x,D);
parameters = get_parameters_jakstat();
lb = parameters.min;
ub = parameters.max;

tolFun  = 1e-10;
tolX    = 1e-12;
maxFunEvals = 2000;
maxIter     = maxFunEvals;
x0 = 0*ones(parameters.number,1);

used_time = cputime;

switch optimizer
    case 'fmincon'
        options = optimset;
        options.MaxFunEvals = maxFunEvals;
        options.MaxIter = maxIter;
        options.Display = 'iter';
        options.Algorithm = 'interior-point';
        options.TolFun = tolFun;
        options.TolX = tolX;
        options.GradObj = 'on';
        options.PrecondBandWidth = Inf;
        [x,fval,exitflag,output] = fmincon(fun,x0,[],[],[],[],parameters.min,parameters.max,[],options);
    case 'scmcr'
        options = struct();
        options.Lb = lb;
        options.Ub = ub;
        options.MaxFunEvals = maxFunEvals;
        options.MaxIter = maxIter;
        [x,fval,meta] = yopt.hoo.scmcr(fun,x0,options);
    case 'Noodles'
        options = struct();
        options.tolGrad = tolX;
        options.tolStep = tolX;
        options.maxObjEvaluations = maxFunEvals;
        options.maxIterations = maxIter;
        results = Noodles(fun,x0,lb,ub,options);
end

used_time = cputime - used_time;

save(fullfile(exdir, [ 'test_jakstat_' optimizer '.mat']));