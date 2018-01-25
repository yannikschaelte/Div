function varargout = ybobyqa(FUN,X0,LB,UB,varargin)
% MatLab wrapper for the mexbobyqa.F routine, which ports Powell's BOBYQA
% routine to MatLab. mexbobyqa.F needs th have been compiled via mex
% before.
%
% Input:
% fun     : objective function to be minimized
% x0      : initial guess for parameters
% lb, ub  : bounds for parameters
% options : struct with options for the algorithm:
%   Rhobeg, Rhoend, MaxFunEvals : see algorithm definition
%
% Output:
% x   : best guess for parameters
% fval: objective function at the solution, generally fval=fun(x)
% exitflag:
%   1 : The function converged to a solution x
%   0 : Number of iterations exceeded options.MaxIter or number of function
%       evaluations exceeded options.MaxFunEvals.
%   -1: The algorithm was terminated inappropriately
% output : struct with meta information:
%   funcCount   : number of function evaluations
%   algorithm   : name of the algorithm
%   t_cpu       : cpu time
%
% Currently has 2nd mode with 0 input arguments returning a function handle
% used by funHandleWrap.m to compute objective function values. This is
% because via mexCallMATLAB only non-anonymous functions (e.g. defined in a
% file) or global function handles can be called, but pesto uses anonymous
% function handles. ATTENTION: This might be not thread-safe.

% check for options
if nargin > 4
    options = varargin{1};
else
    options = struct();
end

X0 = X0(:);
XL = LB(:);
XU = UB(:);

% interpret parameters

N = length(X0);

if (isfield(options,'Npt') && ~isempty(options.Npt))
    NPT		= options.Npt;
else
    NPT  = 2*N+1;
end

if (isfield(options,'Rhobeg') && ~isempty(options.Rhobeg))
    RHOBEG = options.Rhobeg;
else
    RHOBEG = 1e-1;
end

if (isfield(options,'Rhoend') && ~isempty(options.Rhoend))
    RHOEND = options.Rhoend;
else
    RHOEND = 1e-8;
end

if (isfield(options,'MaxFunEvals') && ~isempty(options.MaxFunEvals))
    MAXFUN = options.MaxFunEvals;
else
    MAXFUN = 1000*N;
end

IPRINT = 0; % no output

% track time
starttime = cputime;

% do optimization
[ X,FVAL,FEVALS ] = bobyqa(N,NPT,X0,XL,XU,FUN,RHOBEG,RHOEND,MAXFUN,IPRINT);

% meta information
output.funcCount = FEVALS;
output.algorithm = 'ybobyqa';
output.t_cpu = cputime - starttime;
exitflag = 1;

% return
varargout{1} = X;
varargout{2} = FVAL;
varargout{3} = exitflag;
varargout{4} = output;
end

function [ X,FVAL,FEVALS ] = bobyqa(N,NPT,X,XL,XU,FUN,RHOBEG,RHOEND,MAXFUN,IPRINT)

FEVALS = 0;
CALFUN = @(X) calfun(FUN,X);

% check if NPT acceptable

if NPT < N+2 || NPT > (N+2)*(N+1)/2
    error('ybobyqa: NPT not in required interval');
end

% preliminary adjustments and setup of initial quadratic model

[XBASE,SL,SU,XPT,FVAL,GOPT,HQ,PQ,BMAT,ZMAT] = ...
    prelim(N,NPT,X,XL,XU,RHOBEG);

    function varargout = calfun(fun,args)
        varargout = fun(args);
        FEVALS = FEVALS + 1;
    end

end

function [XBASE,SL,SU,XPT,FVAL,GOPT,HQ,PQ,BMAT,ZMAT]...
    = prelim(N,NPT,NDIM,X,XL,XU,RHOBEG)

% Modify initial X if necessary to avoid conflicts between the bounds and
% the construction of the first quadratic model.
% The lower and upper bounds on moves from the updated X are saved in SL,
% SU.

if any()

% init interpolation variables

XPT = zeros(NPT,N);
BMAT = zeros(NDIM,N);
HQ = zeros((N*(N+1)/2),1);
PQ = zeros(NPT,1);
ZMAT = zeros(NPT,NPT-(N+1));

%

end

function [ BMAT, ZMAT ] = update(N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)

% The arrays BMAT and ZMAT are updated, as required by the new position of
% the interpolation point that has the index KNEW. The vector VLAG has N+NPT
% components, set on entry to the first NPT and last N components of the
% product Hw in equation (4.11) of the Powell (2006) paper on NEWUOA.
% further, BETa is set on entry to the value of the parameter with that
% name, and DENOM is set to the denominator of the updating formula.
% Elements of ZMAT may be treated as zero if their moduli are at most
% ZTEST. The first NDIM elements of W are used for working space.

NPTM = NPT-(N+1);

ZTEST = max(max(abs(ZMAT)));
ZTEST = 1-20*ZTEST;

for J=2:NPTM
    if abs(ZMAT(KNEW,J)) > ZTEST
        TEMP = sqrt(ZMAT(KNEW,1)^2+ZMAT(KNEW,J)^2); 
        TEMPA = ZMAT(KNEW,1)/TEMP;
        TEMPB = ZMAT(KNEW,J)/TEMP;
        for I=1:NPT
            TEMP=TEMPA*ZMAT(I,1)+TEMPB*ZMAT(I,J);
            ZMAT(I,J)=TEMPA*ZMAT(I,J)-TEMPB*ZMAT(I,1);
            ZMAT(I,1)=TEMP;
        end
    end
    ZMAT(KNEW,J)=0;
end

% Put the first NPT components of the KNEW-th column of HLAG into W, and
% calculate the parameters of the updating formula.

for I=1:NPT
    W(I) = ZMAT(KNEW,1)*ZMAT(I,1);
end

ALPHA = W(KNEW);
TAU = VLAG(KNEW);
VLAG(KNEW) = VLAG(KNEW)-ONE;

% Complete the updating of ZMAT.

TEMP = sqrt(DENOM);
TEMPB = ZMAT(KNEW,1)/TEMP;
TEMPA = TAU/TEMP;

for I=1:NPT
    ZMAT(I,1) = TEMPA*ZMAT(I,1)-TEMPB*VLAG(I);
end

% Finally, update the matrix BMAT.

for J=1:N
    JP = NPT+J;
    W(JP) = BMAT(KNEW,J);
    TEMPA = (ALPHA*VLAG(JP)-TAU*W(JP))/DENOM;
    TEMPB = (-BETA*W(JP)-TAU*VLAG(JP))/DENOM;
    for I=1:JP
        BMAT(I,J) = BMAT(I,J)+TEMPA*VLAG(I)+TEMPB*W(I);
        if I > NPT
            BMAT(JP,I-NPT) = BMAT(I,J);
        end
    end
end

end


function 
