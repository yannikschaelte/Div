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
                    
                    % solve subproblem to compute step s
                    this.subproblem.solve();
                    
                    % evaluate value at proposed step
                    fval_new = this.compute_fval(this.state.x + this.subproblem.step);
                    
                    % check step
                    accept_step = this.subproblem.evaluate(fval_new);
                    
                    % update state if demanded and possible
                    if accept_step
                        accept_step = this.update_state(this.state.x + this.subproblem.step);
                    end
                    
                    % adapt subproblem according to whether step accepted
                    this.subproblem.handle_accept_step(accept_step);
                    
                    % check termination criteria
                    this.check_termination();
                end
                
                % output
                this.print_output();
                
            end
            
            results = noodles.NoodleResults(this);
            
        end
        
        function init(this)
            this.exitflag = nan;
            this.state = noodles.NoodleState(this.dim);
            if ~this.update_state(this.init_x)
                this.exitflag = -2;
            end
        end
        
        function fval = compute_fval(this, x)
            fval = this.objfun(x);
            % this.state.feval_count = this.state.feval_count + 1;
        end
        
        function success = update_state(this,x)
            [fval,grad,hess] = this.objfun(x);
            this.state.feval_count = this.state.feval_count + 1;
            
            if ~isfinite(fval) || ~all(isfinite(grad)) || ~all(all(isfinite(hess)))
                success = false;
            else
                fval_old = this.state.fval;
                
                this.state.x    = x;
                this.state.fval = fval;
                this.state.grad = grad;
                this.state.hess = hess;
                this.state.gradnorm = norm(grad, 2);
                this.state.fvaldiff = fval_old - fval;
                success = true;
            end
        end
        
        function check_termination(this)
            this.state.iter_count = this.state.iter_count + 1;
            
            if this.state.gradnorm < this.options.tol_grad
                this.exitflag = 1;
            elseif this.subproblem.stepnorm < this.options.tol_step
                this.exitflag = 2;
            elseif abs(this.state.fvaldiff) < this.options.tol_fvaldiff
                this.exitflag = 3;
            elseif this.state.iter_count > this.options.iter_max ...
                    || this.state.feval_count > this.options.feval_max
                this.exitflag = 0;
            end
        end
        
        function bool_cont = cont(this)
            bool_cont = isnan(this.exitflag);
        end
        
        function print_output(this)
            if this.options.verbosity ~= 0
                if mod(this.state.iter_count,10) == 1
                    fprintf('iter\tfeval\tfval\n');
                end
                fprintf('%d\t%d\t%.6e\n', this.state.iter_count,this.state.feval_count,this.state.fval);
            end
        end
    end
end

