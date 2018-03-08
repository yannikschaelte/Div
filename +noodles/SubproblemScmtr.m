classdef SubproblemScmtr < noodles.NoodleSubproblem
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        
        function this = SubproblemScmtr(options_in)
            if nargin < 1
                options_in = struct();
            end
            
            this.options = noodles.SubproblemScmtr.get_options(options_in);
        end
        
        function solve(this)
            
        end
        
        function accept_step = evaluate(this, objfun)
            
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

