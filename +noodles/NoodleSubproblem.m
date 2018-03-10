classdef (Abstract) NoodleSubproblem < handle
    % The NoodleSubproblem class administers the local search for the next
    % evaluation point. To customize the subproblem solution, derive from
    % this class and implement the various methods used in 
    % NoodleProblem.run_optimization().
    
    properties ( GetAccess = 'public', SetAccess = 'protected' )
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
            % Constructor, called implicitly.
        end
        
        function init(this, noodle_problem)
            % Initialize/clear the subproblem, pull required data from the
            % NoodeProblem.
            %
            % Input:
            % noodle_problem    : super NoodleProblem instance
            
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
            % Update values after state has changed.
            %
            % Input:
            % state   : NoodleState from the super NoodleProblem
            
            this.x    = state.x;
            this.fval = state.fval;
            this.grad = state.grad;
            this.hess = state.hess;
        end
        
    end
    
    methods (Abstract)
        
        % Solve the subproblem and update the step variable so that x+step
        % indicates the predicted best next evaluation point.
        solve(this)
        
        % Determine whether the beforehand computed step should be
        % accepted.
        %
        % Input:
        % fval_new  : objfun(x+step)
        accept_step = evaluate(this, fval_new)
        
        % Update internal state according to whether the last step was
        % accepted or not. Only here relevant variables should be changed,
        % not in evaluate(). The reason for this subdivision is that
        % numerical issues can occur when in the outer routine the
        % derivatives are computed (e.g. they might contain nans or infs).
        handle_accept_step(this, accept_step)
        
    end
    
    methods (Static, Abstract)
        
        % Create an options struct filled with default values, overwrite
        % with user input and check for validity.
        options = get_options(options_in)
        
    end
end

