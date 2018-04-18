classdef TF
% TestFunction container.
    
    properties
        
    end
    
    methods (Static)
        
        function [varargout] = ex1(x)
            % saddle points at (0,0), (0,5), (5,0), global minimum at
            % (5,5). (0,5), (5,0), (5,5) are 2nd-order stationary.
            x1 = x(1);
            x2 = x(2);
            varargout{1} = 1/4*x1^4 + 1/4*x2^4 - 5/3*x1^3 - 5/3*x2^3;
            if nargout > 1
                varargout{2} = [x1^3 - 5*x1^2; x2^3 - 5*x2^2];
                if nargout > 2
                    varargout{3} = [3*x1^2 - 10*x1, 0;
                                                 0, 3*x2^2 - 10*x2];
                end
            end
        end
        
        function [varargout] = ex2(x)
            dim = length(x);
            
            varargout{1} = 0;
            for j = 1:dim
                varargout{1} = varargout{1} + 1/2*j*x(j)^2 - (5*j)*sin(x(j));
            end
            if nargout > 1
                varargout{2} = zeros(dim,1);
                for j = 1:dim
                    varargout{2}(j) = varargout{2}(j) + j*x(j) - (5*j)*cos(x(j));
                end
                if nargout > 2
                    varargout{3} = zeros(dim,dim);
                    for j = 1:dim
                        varargout{3}(j,j) = varargout{3}(j,j) + j + (5*j)*sin(x(j));
                    end
                end
            end
        end
        
        function [varargout] = ex3(x)
            dim = length(x);
            x = x(:);
            
            varargout{1} = (x(1)-2)^2 + 10*sum(x.^2) + 10*(dot(x,x)-1)^2;
            if nargout > 1
                varargout{2} = zeros(dim,1);
                varargout{2}(1) = varargout{2} + 2*(x(1)-2);
                varargout{2} = varargout{2} + 20*x + 40*x*(dot(x,x)-1);
                if nargout > 2
                    varargout{3} = zeros(dim,dim);
                    for j = 1:dim
                        varargout{3}(j,j) = 120*x(j)^2 + 40*(dot(x,x)-1) - 40*x(j)^2 + 20;
                        for k = j+1:dim
                            varargout{3}(j,k) = 80*x(j)*x(k);
                            varargout{3}(k,j) = varargout{3}(j,k);
                        end
                    end
                    varargout{3}(1,1) = varargout{3}(1,1) + 2;
                end
            end
        end
        
    end
end

