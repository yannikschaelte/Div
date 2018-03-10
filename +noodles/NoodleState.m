classdef NoodleState < handle
    
    properties
        x;
        fval;
        grad;
        hess;
        gradnorm;
        fvaldiff;
        iter_count;
        feval_count;
        
        meta
    end
    
    methods
        function this = NoodleState(dim)
            this.x = nan(dim,1);
            this.fval = nan;
            this.grad = nan(dim,1);
            this.hess = nan(dim,dim);
            this.gradnorm = nan;
            this.fvaldiff = inf;
            this.iter_count = 0;
            this.feval_count = 0;
            this.meta = struct();
        end

    end
end

