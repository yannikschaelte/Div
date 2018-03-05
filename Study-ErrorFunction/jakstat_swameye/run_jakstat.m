function [  ] = run_jakstat(shape,outlier,approach)

if nargin < 3
    approach = 'standard';
end

load(['data_jakstat' outlier],'D');
exdir=fileparts(which('run_jakstat.m'));
[parameters,options] = getParametersAndOptions_jakstat(approach);

switch approach
    case 'standard'
        nllh = @(x) nllh_jakstat_gennormal_standard(x,D,shape);
    case 'hierarchical'
        nllh = @(x) nllh_jakstat_hierarchical(x,D,options.sc);
    case 'hierarchical-adjoint'
        nllh = @(x) nllh_jakstat_hierarchical_adjoint(x,D,options.sc);
    case 'hierarchical-offsets'
        nllh = @(x) nllh_jakstat_hierarchical_offsets(x,D,options.sc);
    case 'hierarchical-adjoint-offsets'
        nllh = @(x) nllh_jakstat_hierarchical_adjoint_offsets(x,D,options.sc);
    otherwise
        error('approach not recognized');
end

parameters_res = getMultiStarts(parameters,nllh,options.MS);

save(fullfile(exdir,['results_' num2str(shape) outlier '_' approach '.mat']));

end

