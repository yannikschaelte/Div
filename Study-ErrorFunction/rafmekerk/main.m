%% Preliminaries

close all;
clear;
exdir=fileparts(which('run_rafmekerk.m'));
addpath(fullfile(exdir,'models'));
addpath(fullfile(exdir,'data'));

%% Optimization

shapes = {1,1.1,1.5,1.7,2,3};
nShapes = length(shapes);
for j = 1:nShapes
%     run_rafmekerk(shapes{j},false);
%     run_rafmekerk(shapes{j},true);
end

%% Analyze

cell_outlier = {'','_outlier'};
markers = {'o','+','*','x','s','d','^','v','<','>','p','h','o','+','*','x','s','d','^','v','<','>','p','h','o','+','*','x','s','d','^','v','<','>','p','h'};

% waterfalls
figure('name','waterfalls-rme');
for io = 1:length(cell_outlier)
    outlier = cell_outlier{io};
    for is = 1:nShapes
        shape = shapes{is};
        filename = ['results_' num2str(shape) outlier '_standard.mat'];
        if ~exist(filename,'file')
            continue;
        end
        load(filename,'parameters_res');
        semilogy(parameters_res.MS.logPost,[markers{is} '-'],'DisplayName',[num2str(shape) outlier]);
        hold on;
    end
end
legend('show');

load('data_rafmekerk_noreps.mat','D');

amiOptions.rtol = 1e-10;
amiOptions.atol = 1e-10;
amiOptions.sensi_meth = 'forward';
amiOptions.sensi = 0;

for io = 1:length(cell_outlier)
    outlier = cell_outlier{io};
for is = 1:nShapes
    shape = shapes{is};
    
    filename = ['results_' num2str(shape) '_standard.mat'];
    if ~exist(filename,'file')
        continue;
    end
    load(['results_' num2str(shape) outlier '_standard.mat'],'parameters_res');
    theta = parameters_res.MS.par(:,1);
    fprintf('logPost %s = %.10f\n',num2str(shape),parameters_res.MS.logPost(1));

    amiData.t = D(1).t;
    amiData.Y = D(1).Y(:,:,1);
    amiData.Sigma_Y = nan(size(amiData.Y));
    amiData.condition = [D(1).k(:); shape; gamma(1/shape)];
    amiData = amidata(amiData);

    sol = simulate_rafmekerk_gennormal_standard([], theta, [], amiData, amiOptions);

    for iy = 1:size(D(1).Y,2)
        fig = figure('name',['fit_rme_' num2str(shape) outlier '_' num2str(iy)]);
        plot(D(1).t,sol.y(:,iy));
        hold on;
        plot(D(1).t,D(1).Y(:,iy,1),'.');
        
        ym = D(1).Y(:,iy,1);
        y  = sol.y(:,iy);
        ym_isnan = isnan(ym);
        y_isnan = isnan(y);
        ym(ym_isnan | y_isnan) = [];
        y(ym_isnan | y_isnan) = [];
        ycorr = corr(ym,y);
        fprintf('correlation %s = %.10f\n',num2str(shape),ycorr);
    end
end
end