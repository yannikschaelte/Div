classdef (Abstract) NoodleSubproblem < handle
    
    properties
        dim;
        lb;
        ub;
        
        cur_x;
        cur_grad;
        cur_hess;
        
        meta;
    end
    
    methods
        
        function this = NoodleSubproblem()
            % Constructor
        end
        
        function init(this,noodle_problem)
            if isa(noodle_problem,'NoodleProblem')
                this.dim = noodle_problem.dim;
                this.lb  = noodle_problem.lb;
                this.ub  = noodle_problem.ub;
            else
                error('Invalid input for NoodleOptions.init.');
            end
        end

    end
    methods (Static)
       
        function options = get_options()
            options.tol_step = 1e-10;
            options.tol_grad = 1e-5;
            options.iter_max = Inf;
            options.funeval_max = Inf;
            options.hessian_fcn = 'objective';
            options.verbosity   = 1;
        end
        
    end
end

