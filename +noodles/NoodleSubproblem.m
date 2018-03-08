classdef (Abstract) NoodleSubproblem < handle
    
    properties
        options;
        
        dim;
        
        x;
        grad;
        hess;
        s;
        
        meta;
    end
    
    methods
        
        function this = NoodleSubproblem(options)
            % Constructor
            if nargin < 1
                options = struct();
            end
            this.options = noodles.NoodleSubproblem.get_options(options);
        end
        
        function init(this, noodle_problem)
            if isa(noodle_problem,'NoodleProblem')
                this.dim = noodle_problem.dim;
            else
                error('Invalid input for NoodleOptions.init.');
            end
            
            this.x    = nan(this.dim,1);
            this.grad = nan(this.dim,1);
            this.hess = nan(this.dim,this.dim);
            this.s    = nan(this.dim,1);
            
            this.meta = struct();
        end
        
        function update(this, state)
            this.x    = state.x;
            this.grad = state.grad;
            this.hess = state.hess;
            this.s    = nan(this.dim,1);
        end
        
    end
    
    methods (Abstract)
        
        % solve the subproblem and update the step variable so that x+s
        % indicates the predicted best next point
        solve(this)
        
        % evaluate objfun(x+s) and return whether this step should be
        % accepted
        accept_step = evaluate(this, objfun)
        
        % update internal state according to whether the last step can be
        % accepted or not
        handle_accept_step(this, accept_step)
        
    end
    
    methods (Static)
        
        function options = get_options(options_in)
            options = struct();
            
            % fill from input
            cell_fieldnames = fieldnames(options);
            cell_fieldnames_in = fieldnames(options_in);
            
            for jf = 1:length(cell_fieldnames_in)
                fieldname = cell_fieldnames_in{jf};
                if ~any(strcmp(cell_fieldnames,fieldname))
                    error(['Options field ' fieldname ' does not exist.']);
                end
                options.(fieldname) = options_in.(fieldname);
            end
        end
        
    end
end

