classdef NoodleProblem < handle
    % The NoodleProblem class administers the optimization routine. In the
    % constructor, pass all relevant data, and then call run_optimization()
    % to run the optimization.
    
    properties  ( GetAccess = 'public', SetAccess = 'private' )
        objfun;
        x0;
        dim;
        options;
        subproblem;
        
        start_time;
        state;
        
        exitflag;
        flag_initial;
    end
    
    methods
        function this = NoodleProblem(objfun, x0, options)
            this.objfun     = objfun;
            this.x0         = x0(:);
            this.dim        = size(this.x0,1);
            this.options    = noodles.NoodleOptions(options);
            this.subproblem = this.options.subproblem; % duplicate
        end
        
        function results = run_optimization(this)
            % Run optimization
            % 
            % Output:
            % results    : NoodleResults object
            
            % initialize the optimization problem
            this.init();
            
            % initialize the subproblem
            this.subproblem.init(this);
            
            while this.cont()
                                                
                % output
                this.print();
                
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
                
            end
            
            % print final output
            this.print_final();
            
            results = noodles.NoodleResults(this);
            
        end
        
        function init(this)
            % Initialize/clear all variables
            
            this.start_time = cputime;
            this.exitflag = nan;
            this.state = noodles.NoodleState(this.dim);      
            this.flag_initial = true;
            if length(this.options.lb) == 1
                this.options.lb = this.options.lb * ones(this.dim, 1);
            end
            if length(this.options.ub) == 1
                this.options.ub = this.options.ub * ones(this.dim, 1);
            end
            if ~this.update_state(this.x0)
                this.exitflag = -2;
            end    
        end
        
        function fval = compute_fval(this, x)
            % Wrapper for computing just the function value
            
            fval = this.objfun(x);
            % this.state.feval_count = this.state.feval_count + 1;
        end
        
        function success = update_state(this, x)
            % In the state, update the position, function value and 
            % derivatives
            
            [fval,grad,hess] = this.options.derivative_fcn(this, x);
            this.state.feval_count = this.state.feval_count + 1;
            
            if ~isfinite(fval) || ~all(isfinite(grad)) || ~all(all(isfinite(hess)))
                success = false;
            else
                fval_old = this.state.fval;
                
                this.state.x    = x;
                this.state.fval = fval;
                this.state.grad = grad(:);
                this.state.hess = hess;
                this.state.gradnorm = norm(grad, 2);
                this.state.fvaldiff = fval_old - fval;
                success = true;
            
                this.flag_initial = false;
            end
        end
        
        function check_termination(this)
            % Check all termination conditions and set exitflag accordingly
            
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
            % Return true if optimization should be continued, false if any
            % termination condition has been hit
            
            bool_cont = isnan(this.exitflag);
        end
        
        function print(this)
            % Output function, during optimization.
            
            if this.options.verbosity ~= 0
                if mod(this.state.iter_count,10) == 0
                    fprintf('iter\tfeval\tfval\n');
                end
                fprintf('%d\t%d\t%.6e\n', this.state.iter_count,this.state.feval_count,this.state.fval);
            end
        end
        
        function print_final(this)
            % Final output function, after optimization.
            if this.options.verbosity ~= 0
                fprintf('final value found:\n');
                fprintf('iter\tfeval\tfval\n');
                fprintf('%d\t%d\t%.6e\n', this.state.iter_count,this.state.feval_count,this.state.fval);
            end
        end
              
    end
    
    methods (Static)
        
        function [fval, grad, hess] = objective(problem, x)
            % Take hessian from third output of objective function.
            
            [fval,grad,hess] = problem.objfun(x);
        end
        
        function [fval, grad, hess] = sr1(problem, x)
            % Update hessian via Symmetric Rank 1 update.
            
            [fval,grad] = problem.objfun(x);
            
            if problem.flag_initial
                hess = eye(problem.dim);
            else
                grad_prev = problem.state.grad;
                hess_prev = problem.state.hess;
                step      = problem.subproblem.step;

                y = grad - grad_prev;
                den = (y - hess_prev * step)' * step;
                if abs(den) >= sqrt(eps)
                    num = (y - hess_prev * step) * (y - hess_prev * step)';
                    hess = hess_prev + num / den;
                else
                    hess = hess_prev;
                end
            end
        end
        
        function [fval, grad, hess] = dfp(problem, x)
            % Update hessian via Davidon-Fletcher-Powell formula.
            
            [fval,grad] = problem.objfun(x);
            
            if problem.flag_initial
                hess = eye(problem.dim);
            else
                grad_prev = problem.state.grad;
                hess_prev = problem.state.hess;
                step      = problem.subproblem.step;
                
                y = grad - grad_prev;
                gamma = 1 / (y' * step);
                
                hess = (eye(problem.dim) - gamma * y * step') ...
                    * hess_prev ...
                    * (eye(problem.dim) - gamma * step * y') ...
                    + gamma * (y * y');
            end
        end
        
        function [fval, grad, hess] = bfgs(problem, x)
            % Update hessian via Broyden-Fletcher-Goldfarb-Shanno formula.
            
            [fval,grad] = problem.objfun(x);
            
            if problem.flag_initial
                hess = eye(problem.dim);
            else
                grad_prev = problem.state.grad;
                hess_prev = problem.state.hess;
                step      = problem.subproblem.step;

                y = grad - grad_prev;
                s1 = (y*y') / (y'*step);
                s2 = (hess_prev*(step*step')*hess_prev) / (step'*hess_prev*step);
                
                hess = hess_prev + s1 - s2;
            end
        end
        
        function [fval, grad, hess] = custom_sr1(problem, x)
            % Use exact hessian every problem.dim steps, else sr1 updates
            
            if mod(problem.state.feval_count,problem.dim) == 0
                [fval,grad,hess] = noodles.NoodleProblem.objective(problem,x);
            else
                [fval,grad,hess] = noodles.NoodleProblem.sr1(problem,x);
            end
        end
        
    end
end

