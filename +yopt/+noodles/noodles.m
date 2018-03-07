function results = noodles(fun, x0, options)

% handle options
if nargin < 3
    options = struct();
end
options = NoodleOptions(options);

% initialize problem
problem = NoodleProblem();

% run optimization
results = problem.run_optimization(options);

end