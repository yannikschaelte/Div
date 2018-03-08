classdef NoodleState < handle
    
    properties
        x;
        fval;
        grad;
        hess;
        gradnorm;
        iter_count;
        funeval_count;
    end
    
    methods
        function this = NoodleState(dim)
            this.x = nan(dim,1);
            this.fval = nan;
            this.grad = nan(dim,1);
            this.hess = nan(dim,dim);
            this.gradnorm = nan;
            this.iter_count = 0;
            this.funeval_count = 0;
        end

    end
end

