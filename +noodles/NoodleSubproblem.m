classdef (Abstract) NoodleSubproblem < handle
    
    properties
        options;
        dim;
        lb;
        ub;
        
        cur_x;
        cur_grad;
        cur_hess;
        step;
        proposed_x;
        
        meta;
    end
    
    methods
        
        function this = NoodleSubproblem(options)
            % Constructor
            if nargin < 1
                options = struct();
            end
            this.options = get_options(options);
        end
        
        function init(this,noodle_problem)
            if isa(noodle_problem,'NoodleProblem')
                this.dim = noodle_problem.dim;
                this.lb  = noodle_problem.lb;
                this.ub  = noodle_problem.ub;
            else
                error('Invalid input for NoodleOptions.init.');
            end
            
            this.cur_x = [];
            this.cur_grad = [];
            this.curr_hess = [];
            this.meta = struct();
        end
        
        function update(this, state)
            this.cur_x = state.x;
            this.cur_grad = state.grad;
            this.cur_hess = state.hess;
        end
                  
        function solve(this)
            
        end
        
        function evaluate(this)
            
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

