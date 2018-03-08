function results = noodles(objfun, init_x, options)

% handle options
if nargin < 3
    options = struct();
end

% initialize problem
problem = noodles.NoodleProblem(objfun, init_x, options);

% run optimization
results = problem.run_optimization();

end