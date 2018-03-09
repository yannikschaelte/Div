classdef NoodleOptions < handle
    
    properties
        subproblem = noodles.SubproblemScmtr();
        
        lb  = -inf;
        ub  = inf;
        
        tol_grad        = 1e-6;
        tol_step        = 1e-6;
        tol_fvaldiff    = 1e-6;
        iter_max        = Inf;
        feval_max       = Inf;   
        hessian_fcn     = noodles.NoodleOptions.objective;
        verbosity       = 1;
        
    end
    
    properties (Constant)
        objective   = 'objective';
        sr1         = 'sr1';
    end
    
    methods
        
        function options = NoodleOptions(options_in)
            if isa(options_in,'noodles.NoodleOptions')
                options = options_in;
            elseif isa(options_in,'struct')
                fields_in = fields(options_in);
                for jf = 1:length(fields_in)
                    if ~isprop(options, fields_in{jf})
                        error(['Cannot assign options property ' fields_in{jf} ' to NoodleOptions.']);
                    end
                    options.(fields_in{jf}) = options_in.(fields_in{jf});            
                end
            else
                error('Invalid input for NoodleOptions');
            end
        end
        
    end
end

