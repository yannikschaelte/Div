classdef NoodleResults < handle

    properties ( GetAccess = 'public', SetAccess = 'private' )
        init_x;
        init_fval;
        final_x;
        final_fval;
        final_grad;
        final_hess;
        final_gradnorm;
        iter_count;
        funeval_count;
        cpu_time;
        exitflag;
    end
    
    methods
        
        function results = NoodleResults(problem)
            % Constructor
            if isa(problem,'NoodleProblem')
                results.init_x          = problem.init_x;
                results.init_fval       = problem.init_fval;
                results.final_x         = problem.state.x;
                results.final_fval      = problem.state.fval;
                results.final_grad      = problem.state.grad;
                results.final_hess      = problem.state.hess;
                results.final_gradnorm  = problem.state.gradnorm;
                results.iter_count      = problem.state.iter_count;
                results.funeval_count   = problem.state.funeval_count;
                results.cpu_time        = cputime - problem.start_time;
                results.exitflag        = problem.exitflag;
            else
                error('The NoodleResults constructor needs a NoodleProblem for initialization.');
            end
        end
        
    end
    
end

