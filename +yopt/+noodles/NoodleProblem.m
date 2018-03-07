classdef NoodleProblem < handle
    
    properties
        objfun;
        init_x;
        init_fval;
        dim;
        start_time;
        options;
        lb;
        ub;
        
        state;
        
        accept_step;
        exitflag;
        
    end
    
    methods
        function this = NoodleProblem()
        end
        
        function results = run_optimization(this)
            
            % initialize the optimization problem
            this.initialize();
            
            % define subproblem
            subproblem = NoodleSubproblem(this);
            
            while isnan(this.exitflag)
                
                subproblem.update(
                
                if this.accept_step
                    this.compute_derivatives();
                end
                this.accept_step = false;
                
                % loop over subproblem solving attempts
                while ~this.accept_step && isnan(this.exitflag)
                    
                    % update subproblem using current derivatives
                    subproblem.update(this);
                    
                    % solve subproblem
                    subproblem.solve();
                    
                    % check if proposed step improves objective
                    this.check_proposed_step(subproblem);
                    
                    % check if termination criteria are met
                    this.check_termination();
                end
                
                % compute gradient and hessian
                this.compute_derivatives();
            end
        end
        
                
        function [] = initialize(this)
           this.exitflag = nan;
           this.compute_derivatives();
        end
    end
end

