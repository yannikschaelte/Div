clear;

addpath('../../NOODLES');
addpath('../testmodels');

objfun = @TF.ex1;
x0 = -ones(2,1);

options.subproblem = noodles.SubproblemCr();
options.verbosity = 0;
options.derivative_fcn = @noodles.NoodleProblem.sr1;
result = noodles(objfun,x0,options);
result.final_x
result.final_fval

options = optimoptions(@fminunc);
options.Algorithm = 'trust-region';
options.SpecifyObjectiveGradient = true;
[x,fval] = fminunc(objfun,x0,options)

% observations:
% * scmcr, scmtr better than cr, tr, avoiding the saddle point apparently
% * cr better with sr1 than with real hessian