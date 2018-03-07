function results = noodles(fun, x0, options)

% handle options
if nargin < 3
    options = struct();
end
options = noodles_options();

results = problem.run_optimization(options);

end

function options = noodles_options(options_in)

options = struct();


end