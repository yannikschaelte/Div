classdef NoodleProblem < handle
    
    properties
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
            fval = this.objfun(x);
            % this.state.feval_count = this.state.feval_count + 1;
        end
        
        function success = update_state(this,x)
            switch this.options.hessian_fcn
                case noodles.NoodleOptions.objective
                    [fval,grad,hess] = this.objfun(x);
                case noodles.NoodleOptions.sr1
                    [fval,grad] = this.objfun(x);
                    hess = this.sr1(grad);
                case noodles.NoodleOptions.bfgs
                    [fval,grad] = this.objfun(x);
                    hess = this.bfgs(grad);
            end
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
        
        function print(this)
            if this.options.verbosity ~= 0
                if mod(this.state.iter_count,10) == 0
                    fprintf('iter\tfeval\tfval\n');
                end
                fprintf('%d\t%d\t%.6e\n', this.state.iter_count,this.state.feval_count,this.state.fval);
            end
        end
        
        function print_final(this)
            if this.options.verbosity ~= 0
                fprintf('final value found:\n');
                fprintf('iter\tfeval\tfval\n');
                fprintf('%d\t%d\t%.6e\n', this.state.iter_count,this.state.feval_count,this.state.fval);
            end
        end
       
        function hess = sr1(this, grad)
            if this.flag_initial
                hess = eye(this.dim);
            else
                grad_prev = this.state.grad;
                hess_prev = this.state.hess;
                step      = this.subproblem.step;

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
        
        function hess = bfgs(this, grad)
            if this.flag_initial
                hess = eye(this.dim);
            else
                grad_prev = this.state.grad;
                hess_prev = this.state.hess;
                step      = this.subproblem.step;

                y = grad - grad_prev;
                s1 = (y*y') / (y'*step);
                s2 = (hess_prev*(step*step')*hess_prev) / (step'*hess_prev*step);
                
                hess = hess_prev + s1 - s2;
            end
        end
        
    end
end

