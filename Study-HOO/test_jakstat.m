function [] = test_jakstat(optimizer)

rng(0);

exdir=fileparts(which('test_jakstat.m'));
addpath('../../NOODLES');
addpath('../testmodels/jakstat');
% addpath('../../Noodles');
load('data_jakstat','D');
fun = @(x) nllh_jakstat(x,D);
parameters = get_parameters_jakstat();
lb = parameters.min;
ub = parameters.max;

tolFun  = 1e-10;
tolX    = 1e-12;
maxFunEvals = Inf;
maxIter     = Inf;
x0 = 0*ones(parameters.number,1);

used_time = cputime;

switch optimizer
    case 'fmincon-ip'
        options = optimoptions(@fmincon,'Algorithm','interior-point');
        options.Display = 'iter';
        options.MaxFunctionEvaluations = maxFunEvals;
        options.SpecifyObjectiveGradient = true;
        hessian = @(x,lambda) hessianFcn(fun,x,lambda);
        options.HessianFcn = hessian;
        [x,fval,exitflag,output] = fmincon(fun,x0,[],[],[],[],parameters.min,parameters.max,[],options);
    case 'fmincon-tr'
        options = optimoptions(@fmincon,'Algorithm','trust-region-reflective');
        options.Display = 'iter';
        options.MaxFunctionEvaluations = maxFunEvals;
        options.MaxIterations = maxIter;
        options.HessianFcn = 'objective'; % as 3rd output of fun
        options.SubproblemAlgorithm = 'factorization'; % 'cg' 'factorization'
        options.SpecifyObjectiveGradient = true;
        [x,fval,exitflag,output] = fmincon(fun,x0,[],[],[],[],parameters.min,parameters.max,[],options);
    case 'scmcr'
        options = struct();
        options.Lb = lb;
        options.Ub = ub;
        options.MaxFunEvals = maxFunEvals;
        options.MaxIter = maxIter;
        [x,fval,meta] = yopt.hoo.scmcr(fun,x0,options);
    case 'scmtr_src'
        options = struct();
        options.maxFunEvals = maxFunEvals;
        options.Blub = 1;
        [x,fval,meta] = yopt.hoo.scmtr_src(fun,x0,options);
    case 'scmcr_src'
        options = struct();
        options.maxFunEvals = maxFunEvals;
        [x,fval,meta] = yopt.hoo.scmcr_src(fun,x0,options);
    case 'noodles'
        options.subproblem = noodles.SubproblemCr();
        options.derivative_fcn = @noodles.NoodleProblem.objective;
        result = noodles(fun,x0,options);
    case 'Noodles'
        options = struct();
        options.tolGrad = tolX;
        options.tolStep = tolX;
        options.maxObjEvaluations = maxFunEvals;
        options.maxIterations = maxIter;
        results = Noodles(fun,x0,lb,ub,options);
    otherwise
        error(['Could not recognize optimizer ' optimizer]);
end

used_time = cputime - used_time

save(fullfile(exdir, [ 'test_jakstat_' optimizer '.mat']));

end

function [H] = hessianFcn(fun,x,~)

[~,~,s2nllh] = fun(x);
H = s2nllh;

end
