clear all;
close all;

nExp = 50;
nStarts = 5;

nAlgs = 3;

for k = 1:nAlgs
    res{k} = zeros(nStarts, nExp);
    objfuns{k} = zeros(1, nExp);
    times{k} = zeros(1, nExp);
end

basefile{1} = 'latin-hypercube';
basefile{2} = 'ss-latinHypercube-separatedLHParameters';
basefile{3} = 'ss-latinHypercube-clusteredParameters';

colors{1} = 'b';
colors{2} = 'r';
colors{3} = 'g';

labels{1} = 'old';
labels{2} = 'ss separated';
labels{3} = 'ss clustered';

for j = 1:nExp
    for k = 1:nAlgs
        load(['res/test_js_' basefile{k} '_' num2str(j-1) '_500_5.mat']);
        res{k}(:, j) = parameters_res.MS.logPost(:);
        objfuns{k}(1, j) = nansum(parameters_res.MS.n_objfun) / nStarts;
        times{k}(1, j) = nansum(parameters_res.MS.t_cpu);
    end
end


% best reached function value

for k = 1:nAlgs
    [best_res{k}, index{k}] = sort(res{k}(1, :),'descend');
end

figure;
hold on;
for k = 1:nAlgs
    plot(1:nExp, best_res{k}, [colors{k} '*-'], 'DisplayName', labels{k});
end
title('Best found function value');
legend('location', 'SouthWest');
print('res/best.png');

% number of function evaluations

% objfuns_lh = objfuns_lh(index_lh);
% objfuns_sslh = objfuns_sslh(index_sslh);
for k = 1:nAlgs
    objfuns{k} = sort(objfuns{k}, 'descend');
end

figure;
hold on;
for k = 1:nAlgs
    plot(1:nExp, objfuns{k}, [colors{k} '*-'], 'DisplayName', labels{k});
end
title('Number of function evaluations');
legend;
print('res/fevals.png');

% overall time (approach 2 + few seconds)

for k = 1:nAlgs
    times{k} = sort(times{k}, 'descend');
end

figure;
hold on;
for k = 1:nAlgs
    plot(1:nExp, times{k}, [colors{k} '*-'], 'DisplayName', labels{k});
end
title('Overall time');
legend;
print('res/time.png');

% plot all multistart results

figure;
hold on;
for j = 1:nExp
    for k = 1:nAlgs
        plot(1:nStarts, res{k}(:, j), [colors{k} '*-']);
    end
end
print('res/allstarts.png');