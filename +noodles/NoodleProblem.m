classdef NoodleProblem < handle
    
    properties
        options;
        objfun;
        init_x;
        init_fval;
        dim;
        start_time;
        lb;
        ub;
        
        state;
        subproblem;
        
        exitflag;
        
    end
    
    methods
        function this = NoodleProblem(subproblem,options)
            this.subproblem = subproblem;
            this.options = options;
        end
        
        function results = run_optimization(this)
            
            % initialize the optimization problem
            this.init();
            
            % initialize the subproblem
            this.subproblem.init(this);
            
            while this.cont()              
                
                % update subproblem location
                this.subproblem.update_location(this.state);
                              
                accept_step = false;
                while ~accept_step && this.cont()
                    
                    % solve subproblem
                    this.subproblem.solve();
                    
                    % check proposed step
                    accept_step = this.subproblem.check_proposed_step(@this.compute_fval);
                    
                    if accept_step
                        accept_step = this.update_state(this.subproblem.proposed_x);
                    end
                    this.subproblem.handle_accept_step(accept_step);      
                    
                    % check termination criteria
                    this.check_termination();
                end
                
            end
            
            results = NoodleResults(this);
                
        end
                
        function init(this)
           this.exitflag = nan;
           if ~this.update_state(x)
               this.exitflag = -2;
           end
        end
        
        function fval = compute_fval(x)
            fval = this.objfun(x);          
            % this.state.feval_count = this.state.feval_count + 1;
        end
        
        function success = update_state(this,x)
            [fval,grad,hess] = this.objfun(x);
            this.state.feval_count = this.state.feval_count + 1;
            
            if ~isfinite(fval) || ~all(isfinite(grad)) || ~all(all(isfinite(hess)))
                success = false;       
            else
                this.state.x    = x;
                this.state.fval = fval;
                this.state.grad = grad;
                this.state.hess = hess;
                success = true;
            end
        end
        
        function check_termination(this)
            if this.state.gradnorm < this.options.tol_grad
                this.exitflag = 1;
            elseif this.state.iter_count > this.options.iter_max ...
                    || this.state.feval_count > this.options.feval_max
                this.exitflag = -1;
            end
        end
        
        function bool_cont = cont(this)
            bool_cont = isnan(this.exitflag);
        end
    end
end

