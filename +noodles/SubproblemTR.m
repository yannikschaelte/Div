classdef SubproblemTR < noodles.NoodleSubproblem
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        
        function this = SubproblemTR(options_in)
            if nargin < 1
                options = struct();
            end
            
            this.options = get_options(options_in);
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    
    methods (Static)
        
        function options = get_options(options_in)
            
            
        end
        
    end
end

