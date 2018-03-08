classdef NoodleProblem < handle
    
    properties
        objfun;
        init_x;
        dim;
        options;
        subproblem;
        
        start_time;
        state;
        
        exitflag;      
    end
    
    methods
        function this = NoodleProblem(objfun, init_x, options)
            this.objfun     = objfun;
            this.init_x     = init_x(:);
            this.dim        = size(this.init_x,1);
            this.options    = noodles.NoodleOptions(options);
            this.subproblem = this.options.subproblem; % duplicate
        end
        
        function results = run_optimization(this)
            
            % initialize the optimization problem
            this.init();
            
            % initialize the subproblem
            this.subproblem.init(this);
            
            while this.cont()              
                
                % update subproblem location
                this.subproblem.update(this.state);
                              
                accept_step = false;
                while ~accept_step && this.cont()
                    
                    % solve subproblem
                    this.subproblem.solve();
                    
                    % check proposed step
                    accept_step = this.subproblem.evaluate(@this.compute_fval);
                    
                    % update state if demanded and possible
                    if accept_step
                        accept_step = this.update_state(this.subproblem.proposed_x);
                    end
                    
                    % adapt subproblem according to whether step accepted
                    this.subproblem.handle_accept_step(accept_step);      
                    
                    % check termination criteria
                    this.check_termination();
                end
                
            end
            
            results = NoodleResults(this);
                
        end
                
        function init(this)
           this.exitflag = nan;
           this.state = noodles.NoodleState(this.dim);
           if ~this.update_state(this.init_x)
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
            this.state.iter_count = this.state.iter_count + 1;
            
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

