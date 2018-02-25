function [parameters] = get_parameters_jakstat()

nPar = 17; % maximum number of parameters
minPar = -5*ones(nPar,1);
maxPar = 3*ones(nPar,1);
maxPar(4)  =  6;
maxPar(2)  =  6;
minPar(10) = -6;
minPar(4)  = -3;
minPar(2)  = -3;

parameters.number = nPar;
parameters.min = minPar(1:nPar,1);
parameters.max = maxPar(1:nPar,1);

end
