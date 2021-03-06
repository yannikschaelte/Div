clear all;
close all;

nStarts_nFunEvals = {[5 500], [4 400], [3 333]};

for jGlobal = 1:3
    
    stg = nStarts_nFunEvals{jGlobal};
    
    nExp = 50;
    nStarts = stg(1);
    
    nAlgs = 3;
    
    % isensee d2d polina
    % benchmark collection caro, elba
    % benchmarks sind verf�gbar �ber d2d online, paper, 20 modelle
    % dirichlet approximation of posterior
    % unser kde in hohen dim zu smooth, um vtlg vern�nftig zu approximieren
    
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
            load(['res/test_js_' basefile{k} '_' num2str(j-1) '_' num2str(stg(2)) '_' num2str(stg(1)) '.mat']);
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
    title(['Best found function value, ' num2str(stg(1)) ' runs, ' num2str(stg(1)*stg(2)) ' budget']);
    xlabel('multi-start run');
    ylabel('log-likelihood');
    legend('location', 'SouthWest');
    saveas(gcf, ['res/best_' num2str(stg(2)) '_' num2str(stg(1)) '.png']);
    
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
    xlabel('multi-start run');
    ylabel('# function evaluations');
    legend;
    saveas(gcf, ['res/fevals_' num2str(stg(2)) '_' num2str(stg(1)) '.png']);
    
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
    xlabel('multi-start run');
    ylabel('time');
    legend;
    saveas(gcf, ['res/time_' num2str(stg(2)) '_' num2str(stg(1)) '.png']);
    
    % plot all multistart results
    
    figure;
    hold on;
    for j = 1:nExp
        for k = 1:nAlgs
            plot(1:nStarts, res{k}(:, j), [colors{k} '*-']);
        end
    end
    title('All starts of all multi-start runs');
    xlabel('multi-start run');
    ylabel('# function evaluations');
    saveas(gcf, ['res/allstarts_' num2str(stg(2)) '_' num2str(stg(1)) '.png']);
    
end