addpath('..');

objfun = @(x) sum(x.^2);
init_x = 10*ones(5,1);

result = noodles.noodles(objfun,init_x);