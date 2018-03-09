classdef (Abstract) NoodleSubproblem < handle
    
    properties
        options;
        
        dim;
        lb;
        ub;
        
        x;
        fval;
        grad;
        hess;
        
        step;
        stepnorm;
    end
    
    methods
        
        function this = NoodleSubproblem()
            % Constructor, is called implicitly
        end
        
        function init(this, noodle_problem)    
            this.dim        = noodle_problem.dim;
            this.lb         = noodle_problem.options.lb;
            this.ub         = noodle_problem.options.ub;
            this.fval       = nan;
            this.grad       = nan(this.dim,1);
            this.hess       = nan(this.dim,this.dim);
            this.step       = nan(this.dim,1);
            this.stepnorm   = inf;
        end
        
        function update(this, state)
            this.x    = state.x;
            this.fval = state.fval;
            this.grad = state.grad;
            this.hess = state.hess;
        end
        
    end
    
    methods (Abstract)
        
        % solve the subproblem and update the step variable so that x+s
        % indicates the predicted best next point
        solve(this)
        
        % evaluate objfun(x+s) and return whether this step should be
        % accepted
        accept_step = evaluate(this, fval_new)
        
        % update internal state according to whether the last step can be
        % accepted or not
        handle_accept_step(this, accept_step)
        
    end
    
    methods (Static, Abstract)
        
        options = get_options(options_in)
        
    end
end

