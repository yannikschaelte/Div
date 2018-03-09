classdef SubproblemTr < noodles.NoodlesSubproblem

    properties
        tr_radius
    end
    
    methods
        
        function this = SubproblemTr(options_in)
            if nargin < 1
                options_in = struct();
            end
            
            this.options = noodles.SubproblemTr.get_options(options_in);
        end
        
        function solve(this)
            
        end
        
        function accept_step = evaluate(this, fval_new)
            
        end
        
        function handle_accept_step(this, accept_step)
            
        end

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

