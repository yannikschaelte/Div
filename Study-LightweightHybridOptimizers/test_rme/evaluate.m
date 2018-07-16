nExp = 50;
nStarts = 10;

res_lh = zeros(nStarts, nExp);
res_sslh = zeros(nStarts, nExp);
objfuns_lh = zeros(1, nExp);
objfuns_sslh = zeros(1, nExp);
time_lh = zeros(1, nExp);
time_sslh = zeros(1, nExp);

for j = 1:50
   load(['res/test_rme_latin-hypercube_' num2str(j-1) '_1000_10.mat']);
   res_lh(:, j) = parameters_res.MS.logPost(:);
   objfuns_lh(1, j) = nansum(parameters_res.MS.n_objfun) / nStarts;
   time_lh(1, j) = nansum(parameters_res.MS.t_cpu);
   load(['res/test_rme_ss-latinHypercube-separatedLHParameters_' num2str(j-1) '_1000_10.mat']);
   res_sslh(:, j) = parameters_res.MS.logPost(:);
   objfuns_sslh(1, j) = nansum(parameters_res.MS.n_objfun) / nStarts;
   time_sslh(1, j) = nansum(parameters_res.MS.t_cpu);
end


% best reached function value

[best_res_lh, index_lh] = sort(res_lh(1, :),'descend');
[best_res_sslh, index_sslh] = sort(res_sslh(1, :),'descend');

figure;
plot(1:nExp, best_res_lh, 'r*-');
hold on;
plot(1:nExp, best_res_sslh, 'b*-');
title('Best found function value');
legend('old', 'new', 'location', 'SouthWest');
print('res/best.png');

% number of function evaluations

% objfuns_lh = objfuns_lh(index_lh);
% objfuns_sslh = objfuns_sslh(index_sslh);
objfuns_lh = sort(objfuns_lh, 'descend');
objfuns_sslh = sort(objfuns_sslh, 'descend');

figure;
plot(1:nExp, objfuns_lh, 'r*-');
hold on;
plot(1:nExp, objfuns_sslh, 'b*-');
title('Number of function evaluations');
legend('old', 'new');
print('res/fevals.png');

% overall time (approach 2 + few seconds)

time_lh = sort(time_lh, 'descend');
time_sslh = sort(time_sslh, 'descend');

figure;
plot(1:nExp, time_lh , 'r*-');
hold on;
plot(1:nExp, time_sslh, 'b*-');
title('Overall time');
legend('old', 'new');
print('res/time.png');