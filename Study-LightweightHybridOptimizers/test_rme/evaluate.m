clear all;
close all;

nStarts_nFunEvals = {[10 1000], [2 1000]};

for jGlobal = 1:2

    stg = nStarts_nFunEvals{jGlobal};
    
    nExp = 50;
    nStarts = stg(1);

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
            load(['res/test_rme_' basefile{k} '_' num2str(j-1) '_' num2str(stg(2)) '_' num2str(stg(1)) '.mat']);
            res{k}(:, j) = parameters_res.MS.logPost(:);
            objfuns{k}(1, j) = nansum(parameters_res.MS.n_objfun) / nStarts;
            times{k}(1, j) = parameters_res.time;
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
    print(['res/best_' num2str(stg(2)) '_' num2str(stg(1)) '.png']);

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
    print(['res/fevals_' num2str(stg(2)) '_' num2str(stg(1)) '.png']);
    
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
    print(['res/time_' num2str(stg(2)) '_' num2str(stg(1)) '.png']);
    
    % plot all multistart results

    figure;
    hold on;
    for j = 1:nExp
        for k = 1:nAlgs
            plot(1:nStarts, -log(-res{k}(:, j)), [colors{k} '*-']);
        end
    end
    print(['res/allstarts_' num2str(stg(2)) '_' num2str(stg(1)) '.png']);
    
end