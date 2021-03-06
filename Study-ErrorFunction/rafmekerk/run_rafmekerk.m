function [  ] = run_rafmekerk(shape,outlier,approach)

if nargin < 3
    approach = 'standard';
end

exdir=fileparts(which('run_rafmekerk'));
[parameters,options] = getParametersAndOptions_rafmekerk(approach);

switch approach
    case 'standard'
        if outlier
            load data_rafmekerk_noreps_outlier.mat
            outlier = '_outlier';
        else
            load data_rafmekerk_noreps.mat
            outlier = '';
        end
        nllh = @(x) nllh_rafmekerk_standard(x,D,shape);
    case 'hierarchical'
        load data_rafmekerk.mat
        nllh = @(x) nllh_rafmekerk_hierarchical(x,D,options.sc);
    case 'hierarchical-adjoint'
        load data_rafmekerk.mat
        nllh = @(x) nllh_rafmekerk_hierarchical_adjoint(x,D,options.sc);
    case 'hierarchical-noreps'
        load data_rafmekerk_noreps.mat
        nllh = @(x) nllh_rafmekerk_hierarchical_noreps(x,D,options.sc);
    case 'hierarchical-noreps-adjoint'
        load data_rafmekerk_noreps.mat
        nllh = @(x) nllh_rafmekerk_hierarchical_noreps_adjoint(x,D,options.sc);
    otherwise
        error('approach not recognized');
end

parameters_res = getMultiStarts(parameters,nllh,options.MS);

save(fullfile(exdir,['results_' num2str(shape) outlier '_' approach '.mat']));
end

